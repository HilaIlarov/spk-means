from setuptools import setup, find_packages, Extension
setup(
    name='spkmeans',
    ext_modules=[
        Extension(
            # the qualified name of the extension module to build
            'spkmeans',
            # the files to compile into our module relative to ``setup.py``
            ['spkmeansmodule.c','spkmeans.c']
        )
    ]
)