We typically make binary packages available for new releases. This page contains notes on where to find these packages and how to install them.

For alternative ways to install MDAnalysis, see the page [Install](Install).

## Installing using binary packages (Linux users) ##
As for **MDAnalysis** 0.8.0, binary packages are available for a few Linux distributions. These packages are hosted on the [Open Build Service](http://openbuildservice.org/) instance provided by the [OpenSUSE ](http://www.opensuse.org/) project.

The repository home page is here: https://build.opensuse.org/project/show/home:MDyNMR-FTW:MDAnalysis

Here is a summary on how to use the repository for the supported distribution.

### Debian & Ubuntu ###

On .deb-package-based distributions, like Debian and Ubuntu, you need to add the repository to your sources.list file (most probably located in /etc/apt). The description line for repo follows the following:
```
    deb http://download.opensuse.org/repositories/home:/MDyNMR-FTW:/MDAnalysis/[your_distro_code]  ./
```

Here are the following versions that are supported and their appropriate code:

  * Debian 7.0: Debian\_7.0
  * Ubuntu 12.04: xUbuntu\_12.04
  * Ubuntu 13.10: xUbuntu\_13.10

For example, on a Debian 7 system you should add the following line to your sources.list:
```
    deb http://download.opensuse.org/repositories/home:/MDyNMR-FTW:/MDAnalysis/Debian_7.0 ./
```

### Fedora ###

On Fedora 19, add the following repository to your system:
```
    http://download.opensuse.org/repositories/home:/MDyNMR-FTW:/MDAnalysis/Fedora_19/
```

On Fedora 20, add the following repository to your system:
```
    http://download.opensuse.org/repositories/home:/MDyNMR-FTW:/MDAnalysis/Fedora_20/
```


### OpenSUSE ###

You can install from the software database:
http://software.opensuse.org/package/python-MDAnalysis
