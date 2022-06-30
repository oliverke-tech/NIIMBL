import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="vLab",
    version="0.2.0",
    author="Hua Zheng",
    author_email="zheng.hua1@northeastern.edu",
    description="vLab is a digital twin simulator for N-linked glycosylation cell culture, bireactor, chromatography "
                "simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='LICENSE',
    # url="https://github.com/pypa/sampleproject",
    project_urls={
        # "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=['scripts/raman_test.py'],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    test_suite='src.tests',
    package_data={'': ['data/resources/Network Description.csv']},
    include_package_data=True,
    install_requires=['wheel',
                      'pytest',
                      'setuptools'],
    extras_require={
        'testing': ['pytest'],
    }
)