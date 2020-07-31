# How to get the code

## Send an email to Caoxiang Zhu
  FOCUS is now reposited in GitHub under PrincetonUniversity's GitHub service. For safety reason, it's still private. 
  If you want to have a copy of the code or collaborate together, please send an email to Caoxiang Zhu 
  (czhu[at]pppl.gov or caoxiangzhu[at]gmail.com) with your GitHub account name.
    
  If you are also interested in using PrincetonUniversity's GitHub service, please visit 
  [here](https://www.princeton.edu/researchcomputing/faq/how-do-github-permissions/).
    
## Play with git
 There are millions of tutorials on the internet introducing *how to use git/github*. You can choose anyone you like.
 I personally recommend [a simple guide](http://rogerdudler.github.io/git-guide/index.html). 
 If you are not interested in learning git, here are the basci git commands you might need when using FOCUS.
    
 *First of all, please install git or module load git.*
    
 * **Clone from GitHub**
      
   To get a copy of the code, you need to use the **git clone** command. Here are two ways to use it.
      
   **a. Using https**
  
     In your terminal, please type          
      ```bash
      git clone https://github.com/PrincetonUniversity/FOCUS.git
      ```          
      It may request to verify your GitHub account and password.
      
   **b. Using ssh**
      
      A better way is to use ssh (avoid typing account and password for future using). You can go to this 
      [page](https://help.github.com/articles/connecting-to-github-with-ssh/) setting your ssh-keygen first.
          
      After you have generated a correct keygen, you can type
          
      ```bash
      git clone git@github.com:PrincetonUniversity/FOCUS.git
      ```
      
 * **Pull changes**
    
      You can always get updated by using *pull* command:
      ```bash
      git pull
      ```
      This will update all the branches you have. 
      You can also specify the branch name you want to pull, like `git pull origin develop`. 
      Please see more details talking about [branch](http://nvie.com/posts/a-successful-git-branching-model/).
    
 * **Push changes**
    
      You can also make some modifications and push them to GitHub.
      Here are the procedures.
      
      - First, check modification status by `git status`;
      
      - Then, commit changes to your local server;
      
        You can commit them one by one, like
        
        ```bash
        git add example ; git commit -m "just test"
        ```
        The message behind *-m* is for telling others what changes you made to this file.
        
        You can also commit all the changed files at the same time by
        ```bash
        git commit -am "mutiple commits"
        ```
        
      - Push the changes to GitHub;
      
        If you just want to save on your local server, then `git commit` is enough.
        But if you want others to see your changes, you can push the local commits to GitHub by
        ```bash
        git push origin branch_name
        ```
        The branch name is the branch what you are working on. 
        If you want creat a new branch, type `git branch new` (`git branch -d new` for delete).
        If you want to push the branch to the remote, type `git push origin new` (`git push origin :new for delete`)
        
        More details about branch operations can be seen 
        [here](https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches).
        
## Communications

  The [issues page](https://github.com/PrincetonUniversity/FOCUS/issues) is a good place to discuss and communicate.
  Try it!
  
## *Notes*

  There are no warranties, neither for the code, nor for this page. You are welcome to add/correct contents.
