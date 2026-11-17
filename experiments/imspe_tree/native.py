"""Load the small C++ fitting kernel using SciPy's BLAS/LAPACK ABI."""

from concurrent.futures import ThreadPoolExecutor
import ctypes as ct
from functools import lru_cache
from hashlib import sha256
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from scipy.linalg import cython_blas, cython_lapack


def capsule_address(capsule):
    name = ct.pythonapi.PyCapsule_GetName
    name.argtypes = [ct.py_object]
    name.restype = ct.c_char_p
    pointer = ct.pythonapi.PyCapsule_GetPointer
    pointer.argtypes = [ct.py_object, ct.c_char_p]
    pointer.restype = ct.c_void_p
    return pointer(capsule, name(capsule))


@lru_cache(maxsize=1)
def kernel():
    source = Path(__file__).with_name("lmm_native.cpp")
    signature = sha256(source.read_bytes()).hexdigest()[:16]
    windows = sys.platform == "win32"
    suffix = ".dll" if windows else ".so"
    cache = Path(tempfile.gettempdir()) / "imspe-native" / signature
    cache.mkdir(parents=True, exist_ok=True)
    library = cache / ("lmm" + suffix)
    if not library.exists():
        compiler = os.environ.get("CXX") or shutil.which("g++")
        if compiler is None and windows:
            tool = Path("C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe")
            compiler = str(tool) if tool.exists() else None
        if not compiler:
            raise RuntimeError("Configure CXX for the native solver, or run with --solver torch")
        compiler = shutil.which(compiler) or str(Path(compiler).resolve())
        compiler_env = {**os.environ, "PATH": str(Path(compiler).parent) + os.pathsep + os.environ.get("PATH", "")}
        with tempfile.TemporaryDirectory(dir=cache) as build:
            temporary = Path(build) / library.name
            flags = ["-static-libgcc", "-static-libstdc++"] if windows else ["-fPIC"]
            subprocess.run([compiler, "-std=c++17", "-O3", "-shared", *flags,
                            str(source), "-o", str(temporary)], check=True, env=compiler_env)
            try:
                temporary.replace(library)
            except PermissionError:
                if not library.exists():
                    raise
    dll = ct.CDLL(str(library))
    fit = dll.lmm_fit
    array = np.ctypeslib.ndpointer(dtype=np.float64, flags=("C_CONTIGUOUS", "ALIGNED"))
    fit.argtypes = [array, array, array, ct.c_double,
                    ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_double,
                    ct.c_int, array, array, array] + [ct.c_void_p] * 4
    fit.restype = ct.c_int
    functions = [capsule_address(cython_lapack.__pyx_capi__[name])
                 for name in ("dpotrf", "dpotri", "dpotrs")]
    functions.append(capsule_address(cython_blas.__pyx_capi__["dgemm"]))
    return fit, functions


@lru_cache(maxsize=1)
def pool():
    return ThreadPoolExecutor(max_workers=2, thread_name_prefix="native-ml-em")


def fit_batch(gram, w_k, beta0, ve0, n, rank, tolerance, max_iterations):
    fit, functions = kernel()
    k, q = gram.shape[1], beta0.size

    def one(index):
        beta = np.empty(q)
        effects = np.empty((k, rank))
        diagnostic = np.empty(6)
        status = fit(gram[index], w_k[index], beta0, ve0, q, n, k, rank, tolerance,
                     max_iterations, beta, effects, diagnostic, *functions)
        return (beta, effects, diagnostic) if status == 0 else None

    return list(pool().map(one, range(len(gram)))) if len(gram) > 1 else [one(0)]
