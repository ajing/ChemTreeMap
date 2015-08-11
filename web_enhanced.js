var connect = require('connect');
var http = require('http');
var fs = require('fs');
var qs = require('querystring');
var compression = require('compression')

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

// A help function endwith
function endsWith(str, suffix) {
  return str.indexOf(suffix, str.length - suffix.length) !== -1;
}

function downloadsvg(req, res) {
  if (req.method == "POST"){
      var body = '';
      req.on('data', function (data) {
        body += data;
      });
      req.on('end', function () {
        var POST = qs.parse(body);
        if( endsWith(req.url, "downloadpng") ){
          convertpng(POST.svg, res);
        }
        if( endsWith(req.url, "downloadsvg") ){
          convertsvg(POST.svg, res);
        }
      });
  }
}

var app = connect()
  .use(compression())
  .use(logger)
  .use("/", connect.static(__dirname))
  .use(downloadsvg)
  .listen(8080);

//http.createServer(app)
//    .listen(8080);

//express.use(express.static(__dirname))
//    .post('/Large/download', convertsvg);
