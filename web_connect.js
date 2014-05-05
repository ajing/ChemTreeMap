var connect = require('connect');

function logger(req, res, next) {
  console.log('%s %s', req.method, req.url);
  next();
}

connect()
  .use(logger)

connect.createServer(
    connect.static(__dirname)
).listen(8080);
