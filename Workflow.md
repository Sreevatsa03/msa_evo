# Basic Collaboration Git Workflow:

1. `git pull`

   This updates your main branch with any changes that other people have pushed.

2. `git branch <name>`

   This creates a branch called <name> for you to make useful changes in. Name this branch something simple but memorable, usually your name. For example, I use `git branch sree`.

3. `git checkout <name>`

   This changes which branch you're looking at to the one you just created so that you can make changes without the risk of interfering with other changes to the repo. After running this command, you can start making changes to the files you want to edit.

4. `git add <any files you created or changed>` or `git add .`

   This tells git you want it to commit the changes you've made to the files you've listed. If you choose to do `git add .` instead, git will track changes to all files (use this command for convenience).

5. `git commit` or `git commit -m "commit message"`

   This creates a commit with your changes to the files you've added. When you commit changes, you must include a commit message, a short descrption of the changes made. Running `git commit -m "commit message"` allows you to add your commit message in the command line instead of opening a text editor.

6. `git checkout main`

   This changes which branch you're looking at to main in your local version of the repository -- any changes you've committed to <name> will disappear from the file system (but will reappear when you checkout <name>).

7. `git pull`

   This updates your main branch with any changes that other people have pushed -- we need to do this again because people may have changed things while you were working.

8. `git checkout <name>`

   `git rebase main`

   This appends the changes you've made in <name> to come after whatever other changes people have made in main in a nice line.

9. `git checkout main`

   `git merge <name>`

   If you've done everything correctly, this should be a FAST FORWARD merge. This updates main to include the changes you made in <name>.

10. `git push`

    This pushes the changes you've made to the global repository (in this case, on Github). This will fail if your main branch is not up to date (which occurs when someone pushed changes while you were doing steps 7 or 8). If this fails, you should reset main to before your changes then repeat step 6.

11. `git branch -d <name>`

    This deletes the branch you created to make your changes. We don't need it anymore since those changes are in main. When you want to make more changes, start from step 1. Remember, you can have as many branches as you'd like and name them whatever you want.
