This page outlines how newly [[contributed code|ContributingCode]] is merged into the main repository. 

Note: Please also check the [repository structure](RepositoryStructure) and the [development workflow](DevelopmentWorkflow) for background information.

[git](git) allows us to use a _distributed development model_.

  1. **Anyone** can simply **clone the repository** and work on the clone. (If you have an account on GitHub then just [fork](https://help.github.com/articles/fork-a-repo/) the [MDAnalysis/mdanalysis](https://github.com/MDAnalysis/mdanalysis) repository.)
  1. Once the code in the clone is ready to be integrated in the main development line, **send a "pull" request**. On GitHub, this means to create a [pull request](https://help.github.com/articles/using-pull-requests/). If you have a publicly accessible git repository elsewhere, send a pull request to the [mdnalysis-devel list](http://groups.google.com/group/mdnalysis-devel) so that one of the core developers can review the changes and integrate them.

## Procedure ##
### Repository on GitHub ###
The most convenient approach is to use [pull requests on GitHub](https://help.github.com/articles/using-pull-requests/):

1. Submit the pull request.
2. Any pull request is automatically run through the [[unit tests|UnitTests]] (via travis-ci). If it does not pass, fix and amend you pull request.
3. Developers will ask questions and comment in the pull request. Participate in the discussion.
4. When everything looks good, a developer will merge your pull request into the development version.

### Repository elsewhere ###
If your repository is not on GitHub, the a workflow like the one below will be used:

fullofgrace88 wants to tackle [Issue 2](http://issues.mdanalysis.org/2) so they created the clone [fullofgrace88-mdanalysis-tpr-reader](http://code.google.com/r/fullofgrace88-mdanalysis-tpr-reader/). They will work on the clone. They provide the following information:
```
Please merge my fix for Issue #2 and pull from 

   https://code.google.com/r/fullofgrace88-mdanalysis-tpr-reader/
```

A core developer can create a separate topic branch (here named "Issue/2") in their own local repository for testing anything that comes from the clone:

```
git remote add fullofgrace88 https://code.google.com/r/fullofgrace88-mdanalysis-tpr-reader 
git fetch --all
git checkout -b Issue/2 fullofgrace88/master
```

Once this is done, the developer can switch to the topic branch with
```
git checkout Issue/2
```
(Switching back to the main development line is `git checkout develop`.)

Once a **pull request** comes in, the developer does
```
git checkout Issue/2
git pull
```
and reviews/test the changes. If all goes well, the changes are merged into the development branch with
```
git checkout develop
git merge Issue/2
```
and eventually pushed to the main repository with `git push`.
