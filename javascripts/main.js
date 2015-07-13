//var current_project = location.pathname.split("/")[1];
var current_project = "cgpPindel"

var github_api_url = 'https://api.github.com/repos/cancerit/'.concat(current_project, '/readme?ref=master');

var corsproxy = "http://crossorigin.me/";
var github_url = "https://github.com/cancerit/".concat(current_project, "/raw/master/README.md");
var proxy_url = corsproxy.concat(github_url);

var error_html = "<small>Using fail-over page generation, please notify <a href='mailto:cgp-it@sanger.ac.uk' target='_top'>cgp-it@sanger.ac.uk</a></small><br>";
/*
$.get(final_url)
    .success(function (response) {
	    var converter = new showdown.Converter(),
		html = converter.makeHtml(response);
	    $("#content").html(html);
	})
    .error(function () {
	    var xhr = new XMLHttpRequest();
	    xhr.open('GET', 'https://api.github.com/repos/cancerit/'.concat(current_project, '/readme?ref=master'));
	    xhr.setRequestHeader("Accept", "application/vnd.github.3.html");
	    xhr.setRequestHeader("User-Agent", "CancerIT")
		xhr.send();
	    
	    xhr.onload = function(e){
		$('#content').replaceWith(error_html.concat(xhr.response));
	    }
	});
*/
$.ajax({
	url:github_api_url,
	headers: {
	    "Accept": "application/vnd.github.3.html",
	    "User-Agent": "CancerIT"
	}
    })
    .success(function (response) {
	    $('#content').replaceWith(response); 
	})
    .error(function (r) {
	    $.get(proxy_url)
		.success(function (response) {
			var converter = new showdown.Converter(),
			    html = converter.makeHtml(response);
			$("#content").html(html);
		    });
	});
