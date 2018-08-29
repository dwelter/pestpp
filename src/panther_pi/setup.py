import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="panther_pi",
    version="0.0.1",
    author="Dave Welter",
    author_email="dave@inversemodeler.com",
    description="PANTHER parallel run manger Python interface",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=['panther_pi'],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: MIT License",
        "Operating System :: Windows x64",
    ]
)