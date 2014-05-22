
//RBF
// Y = Wmatrix(Nx1) * Kmatrix(NxN)
// Y's come from the centers, so does the Kmatrix
// so use given points to backsolve weights, then we can
// pump points into the rbf to get surface points

//http://step.polymtl.ca/~rv101/thinplates/

var RBF = function(){
  var centers, ws, ys;
  
  var distance = function(pnt1, pnt2){
    var sum = 0;
    for(var i = 0 ; i < pnt1.length ; i++){
      sum += (Math.pow(pnt1[i] - pnt2[i], 2));
    }
    return Math.sqrt(sum);
  }
  
  //this is going to be a thin-plate spline
  //f(x,y) = a1 + a2x + a3y + SUM(wi * kernel())
  var kernel = function(pnt1, pnt2){
    var r = distance(pnt1, pnt2);
    if(r === 0 ) return 0;
    return Math.pow(r, 2) * Math.log(r);
  }
  
  this.compile = function(cents, y_vals, cb){
   
    setTimeout(function(){
      if(!cents || cents.length === 0){
      // console.error('bad centers array :/');
       cb('bad centers array :/');
       return;
     }
     
     if(!cents.every(function(element){
        return element.length === cents[0].length;
     })) {
        //console.error('centers must have same dimensions :/');
        cb('centers must have same dimensions :/');
        return;
     }
     
     centers = cents;
     ws = [];
     ys = y_vals;
     var matrix = [], matRow = [];
     var P = [], pRow = [];
     for(var i = 0 ; i < centers.length ; i++){
      
      matRow = [];
      pRow = [1];
      for(var k = 0 ; k < centers[i].length; k++){
       pRow.push(centers[i][k]);
      }
      
      for(var j = 0 ; j < centers.length ; j++){
        
        matRow.push(kernel(centers[i], centers[j]));
        
      }
      P.push(pRow);
      matrix.push(matRow.concat(pRow));
      
     }
     
     var pT = $M(P).transpose();
     
     var newRows = pT.elements.map(function(row){
       for(var i = row.length ; i < matrix[0].length; i++){
         row.push(0);
       }
       return row;
     });
     
     for(var i = 0 ; i < newRows.length ; i++){
       matrix.push(newRows[i]);
       ys.push(0);
     }
     
    // console.log("X:");
    // console.dir(matrix);
    // console.log("Y:");
    // console.dir(ys);
     
     ws = this._solve(ys, matrix);
     
     if(!ws){
       cb('rbf failed to compile with given centers./nCenters must be unique :/');
       return;
     }
     
     cb(null, {result: 'success'});
     
    }.bind(this), 0);
  
  };
  
  this._solve = function(b, x){
    //a*x = b
    //a = b * x^-1
    //L = [K p]
    //    [pT 0]
    //L (W | a1 a2 a3) = Y
    var X = $M(x);
    var B = $V(b);
    X = X.inverse();
    return X.multiply(B);
  };
  
  this.getValue = function(pnt){
    var result = 0;
    for(var i = 0 ; i < centers.length; i++){
      result += Number(ws.elements[i]) * kernel(pnt, centers[i]);
    }
    result += Number(ws.elements[centers.length]);
    for(var i = 1 ; i < pnt.length ; i++){
      result += pnt[i] * Number(ws.elements[centers.length+i]);
    }
    return result;
  };
  
  this.getValues = function(pnts, cb){
    setTimeout(function(){
      var resultArr = pnts.map(function(pnt){
          return this.getValue(pnt);
      }, this);
      cb(null, {points: pnts, ys: resultArr});
    }.bind(this), 0);
  };
  
};

//testing

var rbf2D = new RBF();

var pnts2D = [
  [10, 10],
  [20, 10],
  [30, 20],
  [40, 100]
  ];
var testPnts2D = [
  [10, 10],
  [20, 10],
  [30, 20],
  [40, 100]
  ];

rbf2D.compile(pnts2D, [10,20,30,40], function(err, data){
  if(err){
    console.error(err);
    return;
  }
  if(data.result == 'success'){
    console.log('worked!');
    rbf2D.getValues(testPnts2D, function(err, result){
      if(err) {console.error(err); return;}
      
      console.dir(result);
    });
  }
});

var pnts3D = [
  [0, 0, 1],
  [0, 1, 0],
  [1, 0, 0],
  [1, 0, 1]
];
  
var testPnts3D = [
  [0, 0, 1],
  [0, 1, 0],
  [1, 0, 0],
  [1, 0, 1]
];

var rbf3D = new RBF();

rbf3D.compile(pnts3D, [10,10,100,100], function(err, data){
  if(err){
    console.error(err);
    return;
  }
  if(data.result == 'success'){
    console.log('worked!');
    rbf3D.getValues(testPnts3D, function(err, result){
      if(err) {console.error(err); return;}
      
      console.dir(result);
    })
  }
});