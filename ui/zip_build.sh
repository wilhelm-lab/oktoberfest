echo "npm run build;"
npm run build;

cd  ./dist/;
echo "zip -r frontent.zip fonts css js img index.html;"
zip -r frontend.zip fonts css js img index.html;
