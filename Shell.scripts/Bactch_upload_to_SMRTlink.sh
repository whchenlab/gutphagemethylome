API_USER="admin"
API_PASS="admin"
AUTH_TOKEN=$(curl -k -s --user KMLz5g7fbmx8RVFKKdu0NOrJic4a:6NjRXBcFfLZOwHc0Xlidiz4ywcsa -d "grant_type=password&username=$API_USER&password=$API_PASS&scope=sample-setup+run-design+run-qc+data-management+analysis+userinfo+openid+import-dataset" https://${ip}:8243/token | jq -r .access_token)

while read -r infile
do
curl -k -X POST "https://${ip}:8243/SMRTLink/1.0.0/smrt-link/job-manager/jobs/import-dataset" \
    -H "accept: application/json" -H "Content-Type: application/json " -H "Authorization: Bearer $AUTH_TOKEN"\
    -d '{"path":"${FILEPATH}/'${infile}'.subreadset.xml"}'    
done < "file.list.process"
