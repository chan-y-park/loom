# Installation instructions

## Debian 10 Buster (Python 3)

Uses loom branch [stable_v3_3](https://github.com/rda0/loom/tree/stable_v3_3)

System packages:

```bash
apt install apache2
apt install libapache2-mod-wsgi-py3
apt install python3-numpy
apt install python3-matplotlib
apt install python3-scipy
apt install python3-sympy
apt install python3-tk
apt install python3-flask
```

Clone `loom` to your document root:

```bash
git clone https://github.com/rda0/loom.git /var/www/loom
cd /var/www/loom
git checkout stable_v3_3
```

Create a Python 3 virtualenv and install `bokeh`:

```bash
apt install virtualenv
apt install python3-virtualenv
virtualenv -p python3 .venv
source .venv/bin/activate
pip install --no-deps bokeh==0.10.0
deactivate
```

Install SageMath and create a symlink to the `sage` binary:

```bash
mkdir bin
ln -s /path/to/sage bin/sage
```

Apache2 WSGI config:

```apache
<Directory "/var/www/loom">
  <Files "wwwloom.wsgi">
    Require all granted
  </Files>
</Directory>
Alias "/static" "/var/www/loom/static"
<Directory "/var/www/loom/static">
  Require all granted
</Directory>
WSGIDaemonProcess loom user=loom group=loom processes=2 threads=5 display-name=%{GROUP} home=/var/www/loom
WSGIProcessGroup loom
WSGIScriptAlias / /var/www/loom/wwwloom.wsgi
```

Start WSGI website:

```bash
a2enmod wsgi
systemctl restart apache2
```

## Ubuntu 18.04 Bionic (Python 2)

Uses loom branch [stable_v3_2](https://github.com/rda0/loom/tree/stable_v3_2)

System packages:

```bash
apt update
apt upgrade
apt autoremove
apt install python-pip
apt install python-dev
apt install ipython
apt install python-setuptools
apt install pkg-config
apt install libpng-dev
apt install libfreetype6-dev
apt install python-tk
apt install tk-dev
apt install apache2
apt install libapache2-mod-wsgi
apt install sagemath
```

Specific python module versions:

```bash
pip install numpy==1.10.1
pip install matplotlib==1.5.0
pip install scipy==0.15.1
pip install sympy==1.0
pip install bokeh==0.10.0
```

Configure apache2:

```apache
<VirtualHost *:80>
  ServerAdmin webmaster@localhost
  DocumentRoot /var/www
  ErrorLog ${APACHE_LOG_DIR}/error.log
  CustomLog ${APACHE_LOG_DIR}/access.log combined

  WSGIDaemonProcess loom home=/var/www/loom python-path=/var/www/loom
  WSGIScriptAlias /loom /var/www/loom/webmain.wsgi
  <Directory /var/www/loom/>
    WSGIApplicationGroup %{GLOBAL}
    WSGIProcessGroup loom
    WSGIScriptReloading On
    Require all granted
  </Directory>
</VirtualHost>
```

Clone the repo to `/var/www/loom` and checkout branch `stable_v3_2`.

Start WSGI website:

```bash
a2enmod wsgi
systemctl restart apache2
```
