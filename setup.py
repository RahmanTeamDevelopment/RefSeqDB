from setuptools import setup

setup(
    name= 'RefSeqDB',
    version = '0.1.0',
    description = 'A tool for creating RefSeq transcript database',
    url = 'https://github.com/RahmanTeamDevelopment/RefSeqDB',
    author = 'Marton Munz',
    author_email = 'munzmarci@gmail.com',
    license = 'MIT',
    packages=['refseqdb'],
    scripts=['bin/RefSeqDB.py','bin/refseqdb'],
    zip_safe=False
)
