#!/bin/sh


sed -i "s| '\*'; # IE_CORS_ORIGIN| '${CORS_ORIGIN}';|g" /etc/nginx/sites-enabled/default;
sed -i "s|PROXY_PREFIX|${PROXY_PREFIX}|g"               /etc/nginx/sites-enabled/default;


# Make sure the directory for individual app logs exists
mkdir -p /var/log/shiny-server
chown shiny.shiny /var/log/shiny-server

echo "## Galaxy Information ##" >> /etc/R/Renviron
env | grep API | awk '{print "GX_"$1}' >> /etc/R/Renviron 
env | grep URL | awk '{print "GX_"$1}' >> /etc/R/Renviron
env | grep HISTORY_ID | awk '{print "GX_"$1}' >> /etc/R/Renviron
echo "## End Galaxy Information ##" >> /etc/R/Renviron

chmod +x /monitor_traffic.sh
exec /monitor_traffic.sh &
# And nginx in background mode.
exec nginx &

exec shiny-server >> /var/log/shiny-server.log 2>&1
