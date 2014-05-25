var connect = require('connect');
var http = require('http');
var fs = require('fs');
var qs = require('querystring');

function logger(req, res, next) {
  console.log('%s %s', req.method, req.url);
  next();
}

function convertpng(svgsrc, res){

  var samplestream = fs.createWriteStream("sample.png");

  res.setHeader('Content-Disposition', 'attachment;filename=\"new.png\"');
  res.writeHead(200, {'Content-Type': 'application/octet-stream;base64'});
  var convert = require('child_process').spawn("convert", ["svg:", "png:-"]);

  convert.stdout.on('data', function (data) {
    samplestream.write(data);
  });
  convert.on('exit', function(code) {
    samplestream.end();
    res.end();
  });
  convert.stdin.write(svgsrc);
  convert.stdin.end();
}


function convertsvg(svgsrc, res){
  var samplestream = fs.createWriteStream("sample.svg");
  res.setHeader('Content-Disposition', 'attachment;filename=\"new.png\"');
  res.writeHead(200, {'Content-Type': 'application/octet-stream;base64'});
  samplestream.write(svgsrc);
  res.end();
}

function downloadsvg(req, res) {
  if (req.method == "POST"){
      var body = '';
      req.on('data', function (data) {
        body += data;
      });
      req.on('end', function () {
        var POST = qs.parse(body);
        if(req.url == "/Large/downloadpng"){
          convertpng(POST.svg, res);
        }
        if(req.url == "/Large/downloadsvg"){
          convertsvg(POST.svg, res);
        }
      });
  }
}

var app = connect()
  .use(connect.compress())
  .use(logger)
  .use("/", connect.static(__dirname))
  .use(downloadsvg)
  .listen(8080);

//http.createServer(app)
//    .listen(8080);

//express.use(express.static(__dirname))
//    .post('/Large/download', convertsvg);
