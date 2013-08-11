test -z  && . /etc/profile
test -r ~/.alias && . ~/.alias
umask 022
PS1='GRASS 6.2.1 (glc2000):\w > '
PROMPT_COMMAND=/usr/lib/grass/etc/prompt.sh
export PATH="/usr/lib/grass/bin:/usr/lib/grass/scripts:/usr/local/bin:/usr/bin:/bin:/usr/bin/X11:/usr/games"
export HOME="/home/baliola"
