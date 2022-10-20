set -e

docker build -f panss-dti/Dockerfile -t coinstacteam/panss-dti panss-dti/
docker push coinstacteam/panss-dti
docker build -f sans-t1/Dockerfile -t coinstacteam/sans-t1 sans-t1/
docker push coinstacteam/panss-dti
docker build -f panss-t1/Dockerfile -t coinstacteam/panss-t1 panss-t1/
docker push coinstacteam/panss-dti
docker build -f sans-dti/Dockerfile -t coinstacteam/sans-dti sans-dti/
docker push coinstacteam/sans-dti

