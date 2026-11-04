#!/usr/bin/env python3
"""Render the complete Markdown audit to a paginated, searchable CJK PDF."""
import argparse
import html
import json
import os
from pathlib import Path
import subprocess

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate, Frame, KeepTogether, PageTemplate, Paragraph, Spacer, Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[2]
RUNTIME = Path.home() / ".cache/codex-runtimes/codex-primary-runtime/dependencies"
NODE = os.environ.get("NODE", str(RUNTIME / "node/bin/node"))
MARKED = os.environ.get("MARKED_MODULE", str(RUNTIME / "node/node_modules/marked/lib/marked.esm.js"))


def normalize(text):
    return text.translate(str.maketrans({"\u2011": "-", "\u2013": "-", "\u2014": "-", "\u2212": "-", "\u00a0": " "}))


def inline(tokens):
    result = []
    for token in tokens:
        kind = token["type"]
        text = normalize(token.get("text", ""))
        if kind == "strong":
            result.append("<b>" + inline(token["tokens"]) + "</b>")
        elif kind == "em":
            result.append(inline(token["tokens"]))
        elif kind == "codespan":
            result.append('<font color="#25515B">' + html.escape(text) + "</font>")
        elif kind == "link":
            result.append(inline(token.get("tokens", [{"type": "text", "text": text}])))
        elif kind == "br":
            result.append("<br/>")
        elif "tokens" in token:
            result.append(inline(token["tokens"]))
        else:
            result.append(html.escape(text).replace("\n", " "))
    return "".join(result)


class AuditDoc(BaseDocTemplate):
    def afterFlowable(self, flowable):
        if hasattr(flowable, "audit_heading"):
            text, level, anchor = flowable.audit_heading
            self.canv.bookmarkPage(anchor)
            self.canv.addOutlineEntry(text, anchor, level=level, closed=False)


def render(source, target):
    pdfmetrics.registerFont(TTFont("AuditCJK", "/Library/Fonts/Arial Unicode.ttf"))
    pdfmetrics.registerFontFamily("AuditCJK", normal="AuditCJK", bold="AuditCJK", italic="AuditCJK", boldItalic="AuditCJK")
    body = ParagraphStyle("Body", fontName="AuditCJK", fontSize=10, leading=15.2,
                          spaceAfter=7, wordWrap="CJK", alignment=TA_LEFT,
                          allowWidows=0, allowOrphans=0)
    styles = {"body": body}
    styles["h1"] = ParagraphStyle("H1", parent=body, fontSize=18, leading=25, spaceAfter=12, keepWithNext=True)
    styles["h2"] = ParagraphStyle("H2", parent=body, fontSize=13, leading=19, spaceBefore=11, spaceAfter=8,
                                 textColor=colors.HexColor("#164A56"), keepWithNext=True)
    styles["h3"] = ParagraphStyle("H3", parent=body, fontSize=11, leading=17, spaceBefore=9, spaceAfter=6,
                                 textColor=colors.HexColor("#164A56"), keepWithNext=True)
    styles["list"] = ParagraphStyle("List", parent=body, leftIndent=10, firstLineIndent=-7, spaceAfter=5)
    styles["quote"] = ParagraphStyle("Quote", parent=body, leftIndent=10, rightIndent=10,
                                    textColor=colors.HexColor("#234E58"), spaceBefore=5, spaceAfter=8)
    styles["code"] = ParagraphStyle("Code", parent=body, fontSize=9, leading=13,
                                   leftIndent=8, rightIndent=8, backColor=colors.HexColor("#F2F5F6"),
                                   borderPadding=7, spaceBefore=5, spaceAfter=9)
    styles["cell"] = ParagraphStyle("Cell", parent=body, fontSize=9, leading=13, spaceAfter=0)
    styles["cell_header"] = ParagraphStyle("CellHeader", parent=styles["cell"], textColor=colors.white)

    token_script = (
        'import fs from "node:fs"; '
        f'import {{marked}} from {json.dumps(Path(MARKED).as_uri())}; '
        'process.stdout.write(JSON.stringify(marked.lexer(fs.readFileSync(process.argv[1],"utf8"))));'
    )
    tokens = json.loads(subprocess.check_output([NODE, "--input-type=module", "-e", token_script, str(source)], text=True))
    target.parent.mkdir(parents=True, exist_ok=True)
    width, height = A4
    margin = 18 * mm
    body_width = width - 2 * margin
    doc = AuditDoc(str(target), pagesize=A4, leftMargin=margin, rightMargin=margin,
                   topMargin=19 * mm, bottomMargin=19 * mm,
                   title="论文与主实验结果对照 | 2026-10-28", author="RGSS experiment audit")

    def page_frame(canvas, document):
        canvas.saveState()
        canvas.setFont("AuditCJK", 8)
        canvas.setFillColor(colors.HexColor("#4D6066"))
        canvas.drawString(margin, height - 11.5 * mm, "RGSS | 论文与主实验结果对照")
        canvas.drawRightString(width - margin, height - 11.5 * mm, "2026-10-28")
        canvas.setStrokeColor(colors.HexColor("#C8D2D5"))
        canvas.setLineWidth(.4)
        canvas.line(margin, height - 14 * mm, width - margin, height - 14 * mm)
        canvas.setFont("AuditCJK", 7.5)
        canvas.drawString(margin, 10.5 * mm, "基于本轮50次重复；CASE编号按论文对应关系核查")
        canvas.drawRightString(width - margin, 10.5 * mm, f"第 {document.page} 页")
        canvas.restoreState()

    doc.addPageTemplates(PageTemplate(id="Audit", frames=[Frame(margin, 19 * mm,
        body_width, height - 38 * mm, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)], onPage=page_frame))
    heading_count = 0

    def paragraph(token, style="body", prefix=""):
        content = inline(token.get("tokens", [{"type": "text", "text": token.get("text", "")}]))
        return Paragraph(prefix + content, styles[style])

    def flow(blocks):
        nonlocal heading_count
        items = []
        for token in blocks:
            kind = token["type"]
            if kind == "space":
                continue
            if kind == "heading":
                depth = min(token["depth"], 3)
                p = paragraph(token, f"h{depth}")
                heading_count += 1
                p.audit_heading = (normalize(token["text"]), depth - 1, f"heading-{heading_count}")
                items.append(p)
            elif kind in {"paragraph", "text"}:
                items.append(paragraph(token))
            elif kind == "list":
                start = token.get("start", 1)
                for i, item in enumerate(token["items"]):
                    label = f"{start+i}. " if token.get("ordered") else "• "
                    first = True
                    for sub in item["tokens"]:
                        if sub["type"] in {"text", "paragraph"}:
                            items.append(paragraph(sub, "list", html.escape(label) if first else ""))
                            first = False
                        elif sub["type"] != "space":
                            items.extend(flow([sub]))
            elif kind == "blockquote":
                for sub in token["tokens"]:
                    if sub["type"] == "paragraph":
                        items.append(paragraph(sub, "quote"))
                    elif sub["type"] != "space":
                        items.extend(flow([sub]))
            elif kind == "code":
                text = html.escape(normalize(token["text"])).replace("\n", "<br/>")
                items.append(Paragraph(text, styles["code"]))
            elif kind == "table":
                rows = [token["header"]] + token["rows"]
                grid = [[Paragraph(inline(cell["tokens"]), styles["cell_header"] if r == 0 else styles["cell"])
                         for cell in row] for r, row in enumerate(rows)]
                n = len(rows[0])
                if n == 4:
                    widths = [44, 122, 43, body_width - 209]
                elif n == 5 and rows[0][0]["text"] == "模型":
                    widths = [39, 42, 30, (body_width - 111) / 2, (body_width - 111) / 2]
                elif n == 5:
                    widths = [50] + [(body_width - 50) / 4] * 4
                elif n == 3:
                    widths = [60, (body_width - 60) / 2, (body_width - 60) / 2]
                else:
                    widths = [body_width / n] * n
                table = Table(grid, colWidths=widths, repeatRows=1, hAlign="LEFT", spaceBefore=5, spaceAfter=10)
                table.setStyle(TableStyle([
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#285965")),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 6), ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                    ("TOPPADDING", (0, 0), (-1, -1), 6), ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F4F7F8")]),
                    ("LINEBELOW", (0, 0), (-1, -1), .35, colors.HexColor("#D1DCDF")),
                ]))
                items.append(table)
            elif kind == "hr":
                items.append(Spacer(1, 8))
            else:
                raise ValueError(f"Unsupported Markdown block: {kind}")
        return items

    doc.build(flow(tokens))
    print(f"Saved {target}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=ROOT / "result/results_main/PAPER_COMPARISON_2026-10-28.md")
    parser.add_argument("--output", type=Path, default=ROOT / "output/pdf/PAPER_COMPARISON_2026-10-28.pdf")
    options = parser.parse_args()
    render(options.source, options.output)
