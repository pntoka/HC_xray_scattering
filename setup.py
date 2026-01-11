from skbuild import setup

setup(
    name='iobs_ngc',
    version='0.1.2',
    description='Python bindings for iobs_ngc',
    packages=[],  # CMake handles the module
    python_requires='>=3.6',
    install_requires=[],
    cmake_install_dir='.',
    cmake_args=[
        '-DCMAKE_BUILD_TYPE=Release',
    ],
)