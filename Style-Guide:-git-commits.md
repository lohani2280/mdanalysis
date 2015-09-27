The guidelines for committing to the MDAnalysis [[git]] repository are part of the [Style Guide](Style-Guide).


  * Check in sets of changes via [git](git). All changes in one revision should have a common theme. If you implemented two rather different things (say, one bug fix and one new feature) then split your check in into two.
  * _Always_ add a descriptive comment (feel free to be verbose).
    * use a short (<50 characters) subject line that summarizes the change
    * leave a blank line
    * optionally, add additional more verbose descriptions; paragraphs or bullet lists (with `-` or `*`) are good
    * manually break lines at 80 characters
    * manually indent bullet lists
    * See also Tim Pope's [A Note About Git Commit Messages](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html) for rationale.