# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# User specific aliases and functions
alias today='date +%Y%m%d'
alias less='less -R'
alias grep='grep --color=auto'
alias fgrep='fgrep --color=auto'
alias egrep='egrep --color=auto'
alias e='emacs -nw'
alias rc='e ~/.bashrc'
alias sc='unalias -a; source ~/.bashrc'
alias rm='rm -i'
alias clean='rm -vf *~ \#*\# \.*~ 2> /dev/null'
alias ls='ls --color=auto -h'
alias ll='ls -alF'
alias la='ls -la'
alias l='ls -CF'
alias ltr='ls -ltr'
alias ltrl='ls -ltr | less'
alias ltrt='ls -ltr | tail'
alias df='df -h'
alias fs='du -sh'
alias gall='git add -A'
alias gst='git status'
alias gch='git checkout'
alias gri='git rebase -i'
alias grc='git rebase --continue'
alias gcc='git cherry-pick --continue'
alias gbd='git branch -D'
alias gdh='git diff HEAD'
alias gdh1='git diff HEAD~'
alias gname='git diff-tree --no-commit-id --name-only -r'
alias hsha='git log | head -n1 | cut -d" " -f2'
alias starwars='telnet towel.blinkenlights.nl'
