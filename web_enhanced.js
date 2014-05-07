var connect = require('connect');
var http = require('http');
var qs = require('querystring');

function logger(req, res, next) {
  console.log('%s %s', req.method, req.url);
  next();
}

function convertsvg(svgsrc, res){

  //res.writeHead(200, {'Content-Type': 'image/png'});
  var convert = require('child_process').spawn("convert", ["svg:", "png:-"]);

  convert.stdout.on('data', function (data) {
    //console.log(data);
    res.pipe(data);
  });
  convert.on('exit', function(code) {
    res.end();
  });
  convert.stdin.write(svgsrc);
  convert.stdin.end();
}

function downloadsvg(req, res) {
  if (req.method == "POST"){
    var body = '';
    req.on('data', function (data) {
      body += data;
    });
    req.on('end', function () {
      var POST = qs.parse(body);
      // use POST
      console.log(POST.svg);
      convertsvg(POST.svg, res);
    });
  }
}

var app = connect()
  .use(logger)
  .use(connect.compress())
  .use("/", connect.static(__dirname))
 // .listen(8080);

http.createServer(app)
    .use(downloadsvg)
    .listen(8080);
