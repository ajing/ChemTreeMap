var connect = require('connect');

function logger(req, res, next) {
  console.log('%s %s', req.method, req.url);
  next();
}

connect()
  .use(logger);

var app = connect.createServer(
    connect.static(__dirname)
);


  res.writeHead(200, {'Content-Type': 'image/png'});
  var convert = child_proc.spawn("convert", ["svg:", "png:-"]),
      values = (url.parse(req.url, true).query['values'] || ".5,.5")
        .split(",")
        .map(function(v){return parseFloat(v)});

  convert.stdout.on('data', function (data) {
    res.write(data);
  });
  convert.on('exit', function(code) {
    res.end();
  });

    .listen(8080);
