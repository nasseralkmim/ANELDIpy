try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'My Project',
    'author': 'Nasser',
    'url': 'nasseralkmim.github.com',
    'download_url': 'Where to download it.',
    'author_email': 'nasser.alkmim@gmail.com',
    'version': '0.1',
    'install_requires': ['nose', 'numpy', 'scipy', 'matplotlib'],
    'packages': ['eldypy', 'tests'],
    'scripts': [],
    'name': 'eldypy'
}

setup(**config)
