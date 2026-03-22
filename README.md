# `fouvol` - A Fourier proxy for a weakly singular Volterra kernel

---

## 📄 代码解读文档 / Documentation Index

| 文件 | 内容 |
|------|------|
| **[代码逐行解读.md](代码逐行解读.md)** | **全部源代码逐行注释与数学公式对照**（按主函数运行顺序，含行号，共1532行） |
| [第2节公式解读.md](第2节公式解读.md) | 论文第2节所有数学公式的完整推导 |
| [数值方法公式代码解读_Part1.md](数值方法公式代码解读_Part1.md) | 数值方法第一部分：乘积矩形/梯形规则 |
| [数值方法公式代码解读_Part2.md](数值方法公式代码解读_Part2.md) | 数值方法第二部分：Fourier递归求解器 |
| [论文翻译.md](论文翻译.md) | 论文中文翻译 |

> **快速入口**：如果您只需要代码解读，请直接点击上方表格第一行的链接 **[代码逐行解读.md](代码逐行解读.md)**。

---

>Home:
https://github.com/variationalform/fouvol

>Copyright (c) 2020, Simon Shaw
(https://github.com/variationalform, https://www.brunel.ac.uk/people/simon-shaw).

>The moral right of the author has been asserted.

>These codes are free software; you can redistribute them and/or
modify them under the terms of the GNU General Public License Version 3 - the terms of which should accompany this document.

This set of codes can be used to reproduce the results in the paper

>*Approximate Fourier series recursion for problems
involving temporal fractional calculus*

>By: Simon Shaw and John R Whiteman
Brunel University London, 2022

Published in CMAME - **Computer Methods in Applied Mechanics and Engineering**

- <https://doi.org/10.1016/j.cma.2022.115537>
- <https://www.sciencedirect.com/science/article/pii/S0045782522005321?via%3Dihub>
- **_*TO DO:*_** Give volume number/year etc once the paper is in print.

The codes are in python  (with `python 2.7.17` used, as explained in the paper), with a bash script used to manage the batch solution. These scripts live in a folder called `fouvol`. To obtain this folder, `cd` to the directory you want to parent it and type...

```bash
git clone https://github.com/variationalform/fouvol.git
#git clone git@github.com:variationalform/fouvol.git
cd fouvol
./fouvol.py -h
./fouvol.py -v 0 -s 3 -a -0.5 -m 5 -T 10 --T1 0.5 -L 16 --Nt 500 -P
```

This should produce the  help page first, and then in the second command, `png` and `eps` versions of the graphics entitled `varphirep` and  `varphirepextended`. These correspond to the Figures 1 and 2 in Section 1 of the paper. All code corresponding to the creation and storing of `jpg` graphics has been commented out.

Note that LaTeX is used as a `matplotlib` interpreter - if you  don't have LaTeX on the host machine you'll need to go through and comment/alter those lines. 

Next, try this:

```bash
./fouvol.py  -v 0 -s 1 -a -0.5 -T 10 -L 1 --Nt 512 -P
```
This should produce the error in the $N_t = 512$ in the results for the product-Rectangle method with $\alpha = −0.5$ and $T = 10$. Specifically,

```bash
FSerror =  0.0
|u(T)-U| =  0.00120137905322
```
Further,
```bash
./fouvol.py  -v 0 -s 3 -a -0.5 -m 5 -T 10 --T1 0.05 -L 32 --Nt 128
```

should give

```bash
FSerror =  0.123300371378
|u(T)-U| =  0.0253044042786
```
which corresponds to the error in the $N_t = 128$ and $L=32$ Fourier
method with $T = 10$, $\alpha = −0.5$, $T_1 = 0.05$ and a Hermite smoothness of $m = 5$.

The full set of results can be generated using the bash scripts `shortrun.sh` and `bigrun.sh`. These are designed to be run in a linux environment. In this environment execute this

```bash
THEN=`date`
./shortrun.sh | tee shortrun.out && ./compare.sh  | tee -a shortrun.out 
echo $THEN && date
```

to get a stripped-down set of the results in the paper. This should take less than an hour - perhaps around 20 minutes on a fast machine. The full set of results can be obtained by switching these commented lines over in `shortrun.sh`...

```bash
# JARGS="-J 5 16"
# LARGS="-L 3 11"
JARGS="-J 5 7"
LARGS="-L 3 7"
```
Beware though - with `16` and `11` determining the max number of time steps and Fourier components, it will take a while to run...

You can always tidy up with this...

```bash
rm -rf compare_?.sh compare.sh errortable.* *pyc results/ runout.txt timestable.* *.eps *.png *.txt *.out
```
This wont delete important files, but will get rid of everything that can be re-generated with other runs.

## Git management - some notes

Install, with an update first. For example, with Linux Mint:

```bash
sudo apt-get update
sudo apt install git
git config --global user.name "Your Name"
git config --global user.email your@email
git config --global core.editor vi
git config --global color.ui auto
git config --list
git --version
```

As usual, if you make changes, document them and then push as usual. For example,

```bash
git add *
git add .gitignore 
git status
git commit -m 'Initial commit of working code and explanatory README.md'
git push
git pull
git diff --staged
git fetch origin && git diff --name-only main origin/main # or master
# etc etc
```
The push may fail if the clone was made with `https`. If so,

```bash
# Paste ~/.ssh/id_rsa.pub into Git hub web page
# can't use http anymore - as in the clone above, so
ssh -T git@github.com
git remote -v
git remote set-url origin git@github.com:variationalform/fouvol.git
```
REF: <https://stackoverflow.com/questions/17659206/git-push-results-in-authentication-failed>

## Running with `docker`

The short form of the results can also be obtained with `docker`. Do this in a directory you want to share as follows

```bash
newgrp docker # you may not need this
docker pull variationalform/puretime:fouvol
cd shared_directory
docker run -ti --name fouvol -v "$PWD":/home/shared -w \
       /home/fouvol variationalform/puretime:fouvol
# and then execute at the bash prompt any or all of these
source ~/.bashrc

THEN=`date`
./shortrun.sh | tee shortrun.out && ./compare.sh  | tee -a shortrun.out 
echo $THEN && date

rm -rf compare_?.sh compare.sh errortable.* *pyc results/ runout.txt timestable.* *.eps *.png *.txt *.out

```

You `CNTRL-D` to exit. You may also need `chown -R OWNER:GROUP outfile` on any files
`outfile` produced and copied to the  shared folder.

If you exited and then want to re-attach and run again, 

```bash
docker ps --filter "status=exited"
docker start fouvol
docker restart fouvol
docker image ls
docker attach fouvol
source ~/.bashrc
CNTRL-D to exit
```

If you have issues with the docker daemon or group, try these:

```bash
# daemon and group problems
 sudo systemctl start docker
reboot
sudo usermod -a -G docker [user]

grep docker /etc/group
newgrp docker
```

For interest's sake, the docker image was created with the following steps.

```bash
# cd to the directory with the source in it...
# create an image  with this shared directory, pulls if not found locally
docker run -ti --name fouvol -v "$PWD":/home/shared -w /home/fouvol ubuntu
```
This will get you to a `bash` prompt in the image. Here...

```bash
apt-get update
apt-get install python2
apt-get install vim
vi ~/.bashrc 
# and at the bottom type alias python='/usr/bin/python2' and save it
source ~/.bashrc 
# beware that .bashrc wont get read by the login shell.
apt-get install python-pip
pip2 install numpy
pip2 install matplotlib
apt-get install python-tk
pip2 install scipy
```
Now copy over the files from  the share to the image working directory:

```bash
cp ../shared/*.py .
cp ../shared/*.sh .
# create a symlink so 'python' (and not 'python2') works in the scripts   
ln -s /usr/bin/python2 /usr/bin/python 
```

Now, `docker` wont have a `$DISPLAY` variable, and so this fix,
from <https://stackoverflow.com/questions/37604289/tkinter-tclerror-no-display-name-and-no-display-environment-variable>, was introduced into the `python` sources


```bash
import os
# this must come before any other matplotlib imports
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
```

LaTeX is used in `matplotlib` and so we needed this as well,

```bash
apt-get install texlive
apt-get install dvipng texlive-latex-extra texlive-fonts-recommended 

```
this has meant the image is much bigger than it needs to be (is there a workaround?)  

Now put these echo'ed commands at bottom of `~/.bashrc` as a quick reminder of what to do.

```bash
Quickstart Commands:
 rm -rf compare_?.sh compare.sh errortable.* *pyc results/ runout.txt timestable.* *.eps *.png *.txt *.out
 ./shortrun.sh | tee shortrun.out && ./compare.sh  | tee -a shortrun.out
```
and remember to type `source ~/.bashrc` to see them.

The tag was created on <https://hub.docker.com/repository/docker/variationalform/puretime> with the following:

```bash
# get the tag
docker ps --filter "status=exited"
# commit it
docker commit 19d1e9956cbd
# gives
sha256:c8ff6b020b5ca02cb686ba19fbd9613f4ae8cf7e9204a0f4d657fb284707decf
# Then
docker login
docker tag c8ff6b020 variationalform/puretime:fouvol
docker push variationalform/puretime:fouvol
```

Now a downloader can execute and run as explained at the top of this section

Alternatively, once the tar file, made like this,

```bash
docker save c8ff6b020 > fouvol_docker_c8ff.tar
```

is available,

```bash
docker load < fouvol_docker_c8ff.tar
docker run -ti c8ff
```

**_*TO DO:*_** update <https://hub.docker.com/repository/docker/variationalform/puretime> with the paper's DOI

**_*TO DO:*_** make the tarfile available (figshare?)

