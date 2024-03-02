from setuptools import setup, find_packages

setup(
    name='rxnb',  
    version='0.1',  
    packages=find_packages("."),  
    description='Package for calculating reaction barrier',
    long_description=open('README.md').read(),
    author='Li-Cheng Xu', 
    author_email='licheng_xu@zju.edu.cn',  
    url='https://github.com/licheng-xu-echo/RXNBarrier',  
    install_requires=[  
        'rdkit >= 2022.3.3',
        'ase >= 3.22.1',
        'pandas >= 2.0.3',
        'geometric >= 1.0.1'
    ],
    extras_require={
        'conda': [
            'xtb-python',
        ]
    },
    license="MIT",
    python_requires=">=3.8",
    # 其他参数
)