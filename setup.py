from distutils.core import setup, Extension
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    raise ImportError("Requires cython to "
            "be installed before running setup.py (pip install cython)")
try:
    import pysam
except ImportError:
    raise ImportError("Requires pysam to "
            "be installed before running setup.py (pip install pysam)")
try:
    import numpy as np
except ImportError:
    raise ImportError("Requires numpy to "
            "be installed before running setup.py (pip install numpy)")


include_directory = pysam.get_include()
include_directory.extend([np.get_include()])
setup(
    name='filterBamClip',
    version='0.1',
    description='Remove clipped alignments from bam file by some fractions threshold',
    url='',
    author='Douglas Wu',
    author_email='wckdouglas@gmail.com',
    license='MIT',
    packages=['filterBamClip'],
    zip_safe=False,
    scripts = ['bin/filterSoftClip.py'],
    ext_modules = cythonize([Extension(
                            name = 'filterBamClip.bam_splitter',
                            sources = ['filterBamClip/bam_splitter.pyx'],
                            include_dirs = include_directory,
                            )]),
    install_requires=[
          'cython',
          'pysam',
          'numpy'
      ],
    cmdclass = {'build_ext': build_ext}
)
