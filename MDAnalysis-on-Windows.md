MDAnalysis is currently not officially supported on Windows (see Issue [#450](https://github.com/MDAnalysis/mdanalysis/issues/450)). However, some users reported that [MDAnalysis can be installed on Windows 10](https://github.com/MDAnalysis/mdanalysis/issues/1095#issuecomment-263852956) (Anniversary build #14393.447) inside [Bash on Windows](https://msdn.microsoft.com/en-us/commandline/wsl/about) (Ubuntu 14.04 currently; 16.04 works inside current Insider Previews, see below).

## How to install MDAnalysis in *Bash on Ubuntu on Windows*

### Information
Using Windows 10 with the Anniversary Update (Version 1607, Build 14379.5), it's possible to install MDAnalysis inside the  [Windows Subsystem for Linux](https://msdn.microsoft.com/commandline/wsl/about).

With the Anniversary Update, Ubuntu 14.04 is installed. As Python 2.7 never made it into the official repository, one either has to install [Python 2.7 from source](https://www.python.org/downloads/source/) or overcome this issue by installing [Anaconda (Python 2.7)](https://anaconda.org/).
Using a current Insider Preview build of Windows 10 (Build >14901), it's possible to update the Ubuntu installation to 16.04. This comes with a pre-installed Python 2.7.

Having Python 2.7 installed (either the normal library or via Anaconda), one can install MDAnalysis either by ``pip`` or ``conda`` by following the steps described in the [quick start guide](http://www.mdanalysis.org/pages/installation_quick_start/).

### Steps to get MDAnalysis running on Windows 10

1. Follow [this guide](https://msdn.microsoft.com/de-de/commandline/wsl/install_guide) to activate *Bash on Ubuntu on Windows*
2. Having Ubuntu 14.04, install [Anaconda (Python 2.7)](https://anaconda.org/). With Ubuntu 16.04, Python 2.7 is already installed.
3. Follow the [quick start guide](http://www.mdanalysis.org/pages/installation_quick_start/) to install MDAnalysis via ``conda`` or ``pip``.