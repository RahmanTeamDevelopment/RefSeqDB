from setuptools import setup

setup(
    name= 'RefSeqDB',
    version = '0.4.0',
    description = 'A tool for creating RefSeq transcript databases',
    url = 'https://github.com/RahmanTeamDevelopment/RefSeqDB',
    author = 'RahmanTeam',
    author_email = 'rahmanlab@icr.ac.uk',
    license = 'MIT',
    packages=['refseqdb'],
    scripts=['bin/RefSeqDB.py','bin/refseqdb'],
    zip_safe=False
)
