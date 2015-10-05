#!/bin/bash
sudo nginx -c /home/chan/loom/nginx/loom.conf
python test_webagg.py
