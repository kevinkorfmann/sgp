import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="simple_genomic_processor", # Replace with your own username
    version="1.0",
    author="Kevin Korfmann",
    author_email="kevin.korfmann@protonmail.com",
    description="Genomic processing pipeline: From SRA to VCF",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.6',
    install_requires=[
        'pysradb',
        'clint'

        ]
)
