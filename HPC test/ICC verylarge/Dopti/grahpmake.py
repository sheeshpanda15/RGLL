import pyreadr
import matplotlib.pyplot as plt

R_vals = [10, 20, 50]
methods = ['ALL', 'GALL', 'GALLRS']
colors = ['#E63946', '#2A9D8F', '#457B9D']
markers = ['o', 's', '^']

def load_mspe(filepath):
    result = pyreadr.read_r(filepath)
    return result['rec2']

for slope_type in ['fixedslope', 'randomslope']:
    fig, axes = plt.subplots(3, 3, figsize=(13, 11))
    fig.suptitle(f'MSPE — {slope_type}', fontsize=16, fontweight='bold', y=0.98)

    for case_idx in range(1, 10):
        fname = f'C:\\Users\\sheesh\\Desktop\\RGLL\\HPC test\\ICC=0.5\\Iopti\\{slope_type}_case{case_idx}.Rdata'
        df = load_mspe(fname)
        ax = axes[(case_idx - 1) // 3][(case_idx - 1) % 3]

        for method, col, mk in zip(methods, colors, markers):
            vals = df[f'mspe.{method}'].values
            ax.plot(R_vals, vals, color=col, marker=mk, linewidth=2,
                    markersize=7, label=method)

        ax.set_title(f'Case {case_idx}', fontsize=12, fontweight='bold')
        ax.set_xlabel('Number of Groups (R)', fontsize=9)
        ax.set_ylabel('MSPE', fontsize=9)
        ax.set_xticks(R_vals)
        ax.set_xticklabels(['R=10', 'R=20', 'R=50'])
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3, fontsize=11,
               frameon=True, bbox_to_anchor=(0.5, 0.01),
               title='Method', title_fontsize=11)

    plt.tight_layout(rect=[0, 0.05, 1, 0.97])
    plt.savefig(f'C:\\Users\\sheesh\\Desktop\\RGLL\\HPC test\\ICC=0.5\\Iopti\\MSPE_{slope_type}.png', dpi=150, bbox_inches='tight')
    plt.close()