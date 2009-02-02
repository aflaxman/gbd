function plot(container, data) {
  var f = Flotr.draw($(container),
		     data,
		     {mouse:{
		       track: true,
		       color: 'cyan',
		       sensibility: 1, // => distance to show point get's smaller
		       trackDecimals: 2,
		       trackFormatter: function(obj){ return 'x = ' + obj.x +', y = ' + obj.y; }
		       },
		      legend: {
		       position: 'nw',
		       margin: 25,
		      },
		      points: {show: true},
		      lines:  {show: true},
		     }
		    );
};
