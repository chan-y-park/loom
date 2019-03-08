# Installation instructions

## Ubuntu 18.04

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
