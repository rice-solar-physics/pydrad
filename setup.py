from distutils.core import setup

setup(
    name='pydrad',
    license='MIT',
    version='0.1',
    author='Will Barnes',
    url='https://github.com/rice-solar-physics/pydrad',
    package_data={'pydrad': ['configure/templates/*']},
    packages=[
        'pydrad',
        'pydrad.configure',
        'pydrad.parse',
        'pydrad.visualize'
    ],
    author_email='will.t.barnes@gmail.com',
    description='Tools for configuring and parsing HYDRAD simulations'
)
