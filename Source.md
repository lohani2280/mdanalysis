The MDAnalysis source code is stored in a [[git]] repository at https://github.com/MDAnalysis/mdanalysis/.

## Git repository structure ##
The full description of the git repository structure is available [here](RepositoryStructure).

To clone/download/checkout the repository, use one of the following URLs:
1. `git@github.com:MDAnalysis/mdanalysis.git`
2. `https://github.com/MDAnalysis/mdanalysis.git`
3. fork it on GitHub
Option (1) requires that you have a username on GitHub.

### Stable source code ###
Production-ready code is stored in the **[main branch](RepositoryStructure#The_main_branch)**.

You can [browse the master branch](https://github.com/MDAnalysis/mdanalysis).

Please note that main development is done using the **[develop branch](RepositoryStructure#The_develop_branch)**.

### Development source code ###
Development code is stored in the **[develop branch](RepositoryStructure#The_develop_branch)**.

You can [browse the develop branch](https://github.com/MDAnalysis/mdanalysis/tree/develop).

## Source code organisation ##
When you check out the source code you'll get a directory hierarchy that [looks like](https://github.com/MDAnalysis/mdanalysis/tree/develop)
```
mdanalysis/README
           package/
           testsuite/
```
  * [package](https://github.com/MDAnalysis/mdanalysis/tree/develop/package) contains the _MDAnalysis toolkit_ source code.
  * [testsuite](https://github.com/MDAnalysis/mdanalysis/tree/develop/testsuite) contains the code for the [[UnitTests]] and the datafiles

When you [add a new feature](ContributingCode) (in `package`) then you should also add a test case (see [[UnitTests]]) and possibly data files in `testsuite` and commit both together.

(See also [Issue 87](https://github.com/MDAnalysis/mdanalysis/issues/87).)
