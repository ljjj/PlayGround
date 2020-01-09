# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# User specific aliases and functions
alias rc='e ~/.bashrc'
alias sc='source ~/.bashrc'
alias rm='rm -i'
alias delete='mv -t ~/.trash/'
alias clean='delete *~ &> /dev/null; delete #*# &> /dev/null; delete .*~ &> /dev/null'
alias e='emacs -nw'
alias today='date +%Y%m%d'
alias less='less -R'
alias grep='grep --color=auto'
alias ls='ls --color=auto -h'
alias ll='ls -alF'
alias la='ls -la'
alias l='ls -CF'
alias ltr='ls -ltr'
alias ltrt='ls -ltr | tail'
alias df='df -h'
alias fs='du -sh'
alias gall='git add -A'
alias gst='git status'
alias gch='git checkout'
alias gri='git rebase -i'
alias gname='git diff-tree --no-commit-id --name-only -r'
alias hsha='git log | head -n1 | cut -d" " -f2'
alias starwars='telnet towel.blinkenlights.nl'
