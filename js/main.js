function category_to_dist(d,m){
    if('category' in d){
        var c = parseInt(d.category);
        d.category = [];
        for(var i = 0; i < m; i++){
            if(i==c){ d.category.push(1.0); }
            else{     d.category.push(0.0); }
        }
    }
}

//given dist [3.0,4.0,0.0],[h1,h2] => h
function dist_to_color(dist,H){
    var sum = d3.sum(dist),
        d = dist.length,
        a = 0.0;
    for(var i = 0; i < d; i++){
        dist[i] = dist[i]/sum;
        a = a + (i/(d-1)*dist[i]);     //right side 1.0]
    }
    if(H[1]>H[0]){
        return a*(H[1]-H[0])+H[0];
    }else {
        return (1.0-a)*(H[0]-H[1])+H[1];
    }
}

function get_temporal_set(tree){
    var temporal_set = {},
        temporal_list = [];
    function dft(T){
        if(!('children' in T || '_children' in T)){//leaf node
            if('samples' in T.data) {
                for (var j = 0; j < T.data.samples.length; j++) {
                    var temporal = T.data.samples[j].Collection_Date.toString();
                    if (temporal in temporal_set) {
                        temporal_set[temporal] = temporal_set[temporal] + 1;
                    } else {
                        temporal_set[temporal] = 1;
                    }
                }
            }
        }else{
            if('children' in T){
                for(var i = 0; i < T.children.length; i++){
                    dft(T.children[i]);
                }
            }else{
                for(var i = 0; i < T._children.length; i++){
                    dft(T._children[i]);
                }
            }
        }
    }
    dft(tree);
    for(var temporal in temporal_set){
        var date = new Date(temporal);
        var dates = [date.getMonth()+1,date.getDate(),date.getFullYear()];
        temporal_list.push(dates.join('/'));
    }
    return temporal_list.sort(function(a,b){ return (new Date(a)) > (new Date(b)); });
}

function site_center(G){
    var xy = [0.0,0.0],
        n = 0;
    for(var i = 0; i < G.features.length; i++){
        if(G.features[i].geometry.type=='Point') {
            xy[1] = xy[1] + G.features[i].geometry.coordinates[0];
            xy[0] = xy[0] + G.features[i].geometry.coordinates[1];
            n = n + 1.0;
        }
    }
    if(n >= 1.0) {
        xy[0] = xy[0] / n;
        xy[1] = xy[1] / n;
    }
    //console.log(xy);
    return xy;
}

/*take a diatom tree table and calculate the viz tree
  D is a parsed d3.csv data object:
  [{Phylum:p1,Class:c1,Order:o1,Family:f1,Genus:g1,Taxon_name:t1,Short_name:s1,Category:l1},...]
if trim is set, find the straddle point where the tree starts to split and go one back
*/
function phylo_table_to_viz_tree(D,trim){
    var S = {},
        T = {},
        R = [];
    for(var i = 0; i < D.length; i++){
        var r = [];
        Object.keys(D[i]).forEach(function(d){
            r.push(D[i][d]);
            if(D[i][d] in S){ S[D[i][d]]++; }
            else{             S[D[i][d]]=0; }
        });
        R.push(r);
    }

    //make d an object so you can build it our accross several data sets
    function insert(V,row,i,d){
        if(i>=row.length-3){ //leaf node => path=[leaf]
            var t = {'name':row[i],'data':d};
            category_to_dist(t.data,3);
            t.data.descendants=0; //update freq by definition of leaf
            V.children.push(t); //no children on the leaf nodes
        }
        else{ //inner nodes
            if(i<=0){ //root node
                if(!('name' in V)){
                    V.name = row[i];
                    V.data = {'descendants':1};
                    V.children = [];
                }else{
                    V.data.descendants++; //update freq
                }
                insert(V,row,i+1,d);
            }else {
                var v = V.children.map(function(x){ return x.name; }).indexOf(row[i]);
                if(v<0) { //wasn't found among the children
                    var t = {'name': row[i], 'data':{'descendants':1},'children':[]};
                    V.children.push(t); //parent inserts the link to new inner node
                    insert(t,row,i+1,d);
                }else{
                    V.children[v].data.descendants++;  //update freq
                    insert(V.children[v],row,i+1,d);
                }
            }
        }
    }
    //read in the rows of data and update the viz tree structure
    for(var r = 0; r < R.length; r++){
        insert(T,R[r],0,{'category':R[r][R[r].length-1]});
    }
    T.meta = {'n':T.data.descendants};
    if(trim){
        function diverge(T){
            if(!('children' in T)){ //leaf node
                return T; //return just the leaf
            }else{
                if(T.children.length>1){
                    return T;
                }else{
                    return diverge(T.children[0]);
                }
            }
        }
        T = diverge(T);
    }
    return T;
}

/*given a Taxon_name and the sample_table, return all unique sids*/
function phlyo_to_site_ids(st,name){
    var sids = {},
        min = 1.0,
        max = 0.0;
    for(var i = 0; i < st.length; i++){
        if(st[i].Taxon_name==name){
            if(st[i].SID in sids){
                sids[st[i].SID][0] += st[i].RelAbund;
                sids[st[i].SID][1] += 1;
            }else{
                sids[st[i].SID] = [st[i].RelAbund,1];
            }
        }
    }
    for(var sid in sids){
        sids[sid][0] /= 1.0*sids[sid][1]; // average
        if(sids[sid][0]>=max){ max = sids[sid][0]; }
        if(sids[sid][0]<=min){ min = sids[sid][0]; }
    }
    for(var sid in sids){ //normalize [0.0,1.0]
        if(max-min>0.0) {
            sids[sid][0] = (sids[sid][0] - min) / (max - min);
            if (sids[sid][0] < 0.0) {
                sids[sid][0] = 0.0;
            }
            if (sids[sid][0] > 1.0) {
                sids[sid][0] = 1.0;
            }
        }
    }
    return sids;
}

function sample_table_join_viz_tree(S,T,tr,sid,prop){
    var time_set
    function insert(T,name,d){
        if(T.name==name){ //should be a leaf
            if('samples' in T.data){ T.data.samples.push(d); }
            else{ T.data.samples = [d]; }
            T.data.samples.sort(function(a,b){ return a.Collection_Date > b.Collection_Date; }); //ascending by Date(s)
        }else{
            if('children' in T) {
                for (var i = 0; i < T.children.length; i++) {
                    insert(T.children[i],name,d);
                }
            }
        }
    }

    for(var i = 0; i < S.length; i ++){
        S[i].SID = parseInt(S[i].SID);
        S[i].RelAbund = parseFloat(S[i].RelAbund);
        S[i].Collection_Date = new Date(S[i].Collection_Date);
        if(S[i].SID==sid && (S[i].Collection_Date>=tr[0] && S[i].Collection_Date<=tr[1])) {
            insert(T,S[i].Taxon_name,S[i]);
        }
    }

    if(prop){
        function propagation(V,c){
            if(!('children' in V)) { //leaf node
                if('samples' in V.data){
                    V.data.RelAbund = 1e-12;
                    for(var j = 0; j < V.data.samples.length;j++){
                        V.data.RelAbund = V.data.RelAbund+V.data.samples[0].RelAbund;
                    }
                    V.data.RelAbund = V.data.RelAbund/V.data.samples.length;
                } else{
                    V.data.RelAbund = 1e-12;
                }
                return {'RelAbund':V.data.RelAbund,'category':V.data.category};
            }else{         //root or inner nodes
                V.data.RelAbund = 0.0;
                V.data.category = [];
                for(var j = 0; j < c; j++){ V.data.category.push(0.0); }
                for(var i = 0; i < V.children.length; i++){
                    var data = propagation(V.children[i],c);
                    V.data.RelAbund = V.data.RelAbund+data.RelAbund;
                    for(var j = 0; j < c; j++) {
                        V.data.category[j] = V.data.category[j]+data.RelAbund*data.category[j];
                    }
                }
                return {'RelAbund':V.data.RelAbund,'category':V.data.category};
            }
        }
        propagation(T,3); //3 categories => array of 3 => [1.0,0.0,0.0]
    }
    //insert some meta data on a sorted time set to graph.
    return T;
}

function phlyo_tree_map_graph(phylo_data_url,sample_data_url,phylo_id,
                              site_geojson_url,map_id,mapbox) {
    var colors = {0:'rgba(90,180,90,0.5)',1:'rgba(216,179,101,0.5)',2:'rgba(245,245,245,0.5)'}; //can be anything you want
    var color_range = [250,50]; //out of [0,360]

    var attach_id = document.getElementById('phis_viz');
    var margin = {top:   attach_id.clientHeight/45,
                  bottom:attach_id.clientHeight/45,
                  right: attach_id.clientWidth/45,
                  left:  attach_id.clientWidth/45};
    var tree_width_factor = 1.3;
    var width = attach_id.clientWidth/tree_width_factor - (margin.left + margin.right);
    var height = attach_id.clientHeight - (margin.top + margin.bottom);

    var i = 0,
        duration = 300,
        root = {'x0':height/2,'y0':0},
        temporal_list = [],
        phylo_table = [],
        sample_table = [];
    var diagonal = d3.svg.diagonal()
        .projection(function (d) {
            return [d.y, d.x];
        });

    var svg = d3.select('#'+phylo_id).append("svg")
        .attr("width", width + margin.right + margin.left)
        .attr("height", height + margin.top + margin.bottom)
        /*.call(d3.behavior.zoom().on("zoom",function(){
            var xy = d3.event.translate;
            xy[0] = xy[0] + 2*margin.left;
            xy[1] = xy[1] + 2*margin.top;
            svg.attr("transform","translate("+xy[0]+","+xy[1]+")"+" scale("+d3.event.scale+")")
        }))*/
        .append("g")
        .attr("transform", "translate(" + 2*margin.left + "," + margin.top + ")");

    /*
        dt is the diaton tree structure, st is the site table data structure to insert into the dt
        tr is the time_range object to query in the st and the sid is the sid to query in the st
     */
    function refresh_tree(pt,st,tr,sid){
        var x0 = root.x0;
        var y0  = root.y0;
        root = phylo_table_to_viz_tree(pt,true);
        root = sample_table_join_viz_tree(st,root,tr,sid,true); //could select right here...
        root.x0 = x0;
        root.y0 = y0;
        function collapse(d) {
            if (d.children) {
                d._children = d.children;
                d._children.forEach(collapse);
                d.children = null;
            }
        }
        //root.children.forEach(collapse); //collapse the children
        update(root);
    }

    function refresh_temporal_controls(pt,st,tl,sid){
        last_time_slice = null;
        selected_time_slice = null;
        var span = $('<span class="temporal_controls"></span>').appendTo('#temporal_bar');
        var selected_style   = {'color':'#fff','border-color':'#fff','background-color':'rgba(30,30,30,0.5)'};
        var deselected_style = {'color':'#000','border-color':'#000','background-color':'rgba(200,200,200,0.5)'};
        for(var i = 0; i < tl.length; i++){
            $('<button id="temporal_button_'+i+'">'+tl[i]+'</button>').appendTo(span);
            d3.select('#temporal_button_'+i).datum(new Date(tl[i]))
                .on("click", function(d){
                    var t_sid = sid;
                    selected_time_slice = d3.select(this);
                    if(last_time_slice!=null){ //set to null when you want to toggle off
                        if(last_time_slice[0][0]===selected_time_slice[0][0]){ //toggle off
                            console.log('last time slice is toggling off')
                            last_time_slice.style(deselected_style);
                            last_time_slice = null;
                            t_sid = -1;
                        }else{ //new time slice was selected
                            console.log('new time slice was selected, last time slice is toggling off')
                            last_time_slice.style(deselected_style);
                            selected_time_slice.style(selected_style);
                            last_time_slice = selected_time_slice;
                        }
                    }else{//last_time_slice === null => no time slice is on
                        console.log('no time slice is selected');
                        selected_time_slice.style(selected_style);
                        last_time_slice = selected_time_slice;
                    }
                    refresh_tree(pt,st,[d,d],t_sid);
                });
        }

        //set the first temporal button as selected
        last_time_slice = d3.select('#temporal_button_'+0);
        last_time_slice.style(selected_style);
        refresh_tree(pt,st,[new Date(tl[0]),new Date(tl[0])],sid);
    }

    var time_range = [new Date('01/01/2012'), new Date('01/01/2018')],
        tree_width_factor = 1.3,
        geojsonFeature = null,
        lastsitesLayer = null,
        selected_geo_json = {'last':null},
        selected_time_slice = null,
        last_time_slice = null,
        selected_phis_node = null,
        last_phis_color = null;
        mymap = null;
    $.getJSON(site_geojson_url,function(geo_data) {
        d3.csv(phylo_data_url, function (phylo_error, phylo_data) {
            d3.csv(sample_data_url, function (sample_error, sample_data) {
                if (phylo_error) throw phylo_error;
                if (sample_error) throw sample_error;
                phylo_table = phylo_data;
                sample_table = sample_data;

                //map time-----------------------------------------
                geojsonFeature = geo_data; //load your geojson data
                var attach_id = document.getElementById('phis_viz');
                var width = 3 * (attach_id.clientWidth / 2.0);
                var height = attach_id.clientHeight;
                $(map_id).css({'height': height, 'width': width});
                // console.log(geojsonFeature);

                mymap = new L.Map(map_id);
                mymap.setView(site_center(geojsonFeature), 9)
                mymap.doubleClickZoom.disable();
                //basemap loading via OSM or mapbox
                if (mapbox) {
                    L.tileLayer('https://api.tiles.mapbox.com/v4/{id}/{z}/{x}/{y}.png?access_token=pk.eyJ1IjoibWFwYm94IiwiYSI6ImNpejY4NXVycTA2emYycXBndHRqcmZ3N3gifQ.rJcFIG214AriISLbB6B5aw', {
                        maxZoom: 16,
                        minZoom: 9,
                        attribution: 'Map data &copy; <a href="http://openstreetmap.org">OpenStreetMap</a> contributors, ' +
                        '<a href="http://creativecommons.org/licenses/by-sa/2.0/">CC-BY-SA</a>, ' +
                        'Imagery © <a href="http://mapbox.com">Mapbox</a>',
                        id: 'mapbox.streets'
                    }).addTo(mymap);
                } else {
                    // create the tile layer with correct attribution
                    var osmUrl = 'http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png';
                    //var osmUrl = 'https://tiles.wmflabs.org/bw-mapnik/${z}/${x}/${y}.png';
                    var osmAttrib = 'Map data © <a href="https://openstreetmap.org">OpenStreetMap</a> contributors';
                    var osm = new L.TileLayer(osmUrl, {minZoom: 9, maxZoom: 16, attribution: osmAttrib});
                    mymap.addLayer(osm);
                }

                var selected_options = {
                    radius: 12,
                    fillColor: "rgb(250,250,250)",
                    color: "rgb(0,150,250)",
                    opacity: 0.8,
                    fillOpacity: 0.25
                };

                var unselected_options = {
                    radius: 6,
                    fillColor: "rgb(25,25,25)",
                    color: "rgb(200,200,200)",
                    opacity: 0.6,
                    fillOpacity: 0.75
                };

                var hover_options = {
                    radius: 8,
                    fillColor: "rgb(75,75,75)",
                    color: "rgb(0,150,250)",
                    opacity: 0.8,
                    fillOpacity: 0.75
                };

                function onEachFeature(feature, layer) {
                    layer.on('click',function(e){
                        var site = e.target.feature.properties;
                        if(selected_geo_json['last']==null){
                            //console.log('first selection')
                            e.target.setStyle(selected_options);
                        }else{
                            var last = selected_geo_json['last'].target.feature.properties
                            if(site.SID!=last.SID){
                                //console.log('new selection')
                                selected_geo_json['last'].target.setStyle(unselected_options);
                                e.target.setStyle(selected_options);
                            }else{ //already have it selected
                                //console.log('already was selected')
                                if(e.target.options.fillColor==selected_options.fillColor) {
                                    e.target.setStyle(unselected_options);
                                    site = {'SID':-1,'Station_Name':''};
                                }else{
                                    e.target.setStyle(selected_options);
                                }
                            }
                        }
                        selected_geo_json['last'] = e;
                        refresh_tree(phylo_table,sample_table,time_range,site.SID);
                        if(site.SID>=0) {
                            var nav = $('#nav_site').text('SID=' + site.SID + ' : ' + site.Station_Name.toUpperCase());
                        }else{
                            var nav = $('#nav_site').text('phis_viz');
                        }
                        //setup the temporal bar, make into seperate function...
                        temporal_list = get_temporal_set(root);
                        //console.log(temporal_list);
                        $('#temporal_bar').html('');
                        if(temporal_list.length>0) { //make a set of date buttons that will invoke samples
                            refresh_temporal_controls(phylo_table,sample_table,temporal_list,site.SID);
                        }
                    });

                    layer.on('mouseover',function(e){
                        var site = e.target.feature.properties;
                        var last = {'SID':-1};
                        var nav_text = $('#nav_site').text();
                        //console.log(nav_text)
                        if(selected_geo_json['last']!=null) {
                            last = selected_geo_json['last'].target.feature.properties;
                        }
                        if(site.SID!=last.SID) {
                            e.target.setStyle(hover_options);
                            var nav = $('#nav_site').text(nav_text+' > '+site.Station_Name.toUpperCase());
                        }
                    });

                    layer.on('mouseout',function(e){
                        var site = e.target.feature.properties;
                        var last = {'SID':-1};
                        var nav_text = $('#nav_site').text();
                        if(selected_geo_json['last']!=null) {
                            last = selected_geo_json['last'].target.feature.properties;
                        }
                        if(site.SID!=last.SID){
                            e.target.setStyle(unselected_options);
                            var nav_text = $('#nav_site').text($('#nav_site').text().split(' > ')[0]);
                        }
                    });
                }
                L.geoJSON(geojsonFeature,{
                    pointToLayer: function (feature, latlng) {
                        return L.circleMarker(latlng, unselected_options);
                    },
                    onEachFeature: onEachFeature
                }).addTo(mymap);

                lastsitesLayer = L.geoJSON({type: "FeatureCollection",features:[]},{
                    pointToLayer: function (feature, latlng) {
                        return L.circleMarker(latlng, unselected_options);
                    },
                    onEachFeature: onEachFeature
                }).addTo(mymap);

                mymap.on('zoomend', function() {
                    var zoom_value = mymap.getZoom();
                    //console.log(Object.keys(L));
                });
                refresh_tree(phylo_table,sample_table,time_range,-1); //init is a 0.0 relabund blank tree
            });
        });
    });

    $('#map_toggle_button').click(function(){
        var w = parseInt($('#map_side').css('width').split('px')[0]);
        var v = document.getElementById('phis_viz').clientWidth;
        if(w<=0.1){
            console.log('map is open...');
            tree_width_factor = 1.3;
            $('#map_side').css('width','49%');
        }
        if(w/v > 0.45 && w/v < 0.55){
            console.log('map is closed...');
            tree_width_factor = 1.0;
            $('#map_side').css('width','0%');
        }
        update(root);
        //this might have to redraw the selected node as well...
    })

    function update(source) {
        var margin = {top:   attach_id.clientHeight/45,
                      bottom:attach_id.clientHeight/45,
                      right: attach_id.clientWidth/45,
                      left:  attach_id.clientWidth/45};
        var w = parseInt($('#map_side').css('width').split('px')[0]);
        //console.log(w);
        if(w<=1){
            var width = attach_id.clientWidth - (margin.left + margin.right);
        }else {
            var width = attach_id.clientWidth/tree_width_factor - (margin.left + margin.right);
        }
        var height = attach_id.clientHeight - (margin.top + margin.bottom);
        var tree = d3.layout.tree().size([height,width+margin.left+margin.right]);
        d3.select(phylo_id).select('svg')
            .attr("width", width + margin.right + margin.left)
            .attr("height", height + margin.top + margin.bottom);

        // Compute the new tree layout.
        var nodes = tree.nodes(root).reverse(),
            links = tree.links(nodes);
        // Normalize for fixed-depth.
        nodes.forEach(function (d) {
            d.y = d.depth * width/8;
        });
        // Update the nodes…
        var node = svg.selectAll("g.node")
            .data(nodes, function (d) {
                return d.id || (d.id = ++i);
            });
        // Enter any new nodes at the parent's previous position.
        var nodeEnter = node.enter().append("g")
            .attr("class", "node")
            .attr("transform", function(d){
                return "translate(" + source.y0 + "," + source.x0 + ")";
            })
            .attr('name',function(d){ return d.name });
        nodeEnter.append("circle")
            .attr("r", 1e-6)
            .on('mouseover',mouse_over)
            .on('mouseout',mouse_out)
            .on("click", left_click)        //toggle children off and on
        //.on('contextmenu', d3.contextMenu(menu));  //context menu or add node?

        nodeEnter.append("text")
            .attr("x", function (d) {
                return d.children || d._children ? 0 : 10;
            })
            .attr("y", function (d) {
                return d.children || d._children ? -1*(d.data.RelAbund*50+10) : 0;
            })
            .attr("dy", ".35em")
            .style('font-size','0.75vw')
            .attr("text-anchor", function (d) {
                return d.children || d._children ? "end" : "start";
            })
            .text(function (d) {
                var stub = d.name.split(' ');
                if(stub.length>1) { stub = stub.splice(1).join(' '); }
                else{ stub = stub[0]; }
                return d.children ? '' : stub;
            })
            .style("fill-opacity", 1e-6);

        // Transition nodes to their new position.
        var nodeUpdate = node.transition()
            .duration(duration)
            .attr("transform", function (d) {
                return "translate(" + d.y + "," + d.x + ")";
            });

        nodeUpdate.select("circle")
            .attr("r", function(d){
                return d._children ? 50*d.data.RelAbund+8:40*d.data.RelAbund+6;
            })
            .style("fill", function(d) {
                var c = 'hsla(' + dist_to_color(d.data.category, color_range) + ',100%,50%,0.75)';
                if (d.data.RelAbund < 1e-9){
                    c = 'rgba(255,255,255,0.5)';
                }
                if(selected_phis_node!=null){
                    var selected = d3.select(d3.select(selected_phis_node).node().parentNode);
                    var selected_name = selected.attr('name');
                    if(d.name==selected_name) {
                        c = 'hsla(340,100%,50%,0.6)';
                    }
                }
                return d._children ? "rgba(0,0,0,0.5)":c;
            })

        nodeUpdate.select("text")
            .attr("x", function (d) {
                return d.children || d._children ? 0 : 10;
            })
            .attr("y", function (d) {
                return d.children || d._children ? -1*(d.data.RelAbund*50+10) : 0;
            })
            .attr("dy", ".35em")
            .attr("text-anchor", function (d) {
                return d.children || d._children ? "end" : "start";
            })
            .style("fill-opacity", 1);

        // Transition exiting nodes to the parent's new position.
        var nodeExit = node.exit().transition()
            .duration(duration)
            .attr("transform", function (d) {
                return "translate(" + source.y + "," + source.x + ")";
            })
            .remove();

        nodeExit.select("circle")
            .attr("r", 1e-6);
        nodeExit.select("text")
            .attr("x", function (d) {
                return d.children || d._children ? 0 : 10;
            })
            .attr("y", function (d) {
                return d.children || d._children ? -1*(d.data.RelAbund*50+10) : 0;
            })
            .attr("dy", ".35em")
            .attr("text-anchor", function (d) {
                return d.children || d._children ? "end" : "start";
            })
            .style("fill-opacity", 1e-6);

        // Update the links…
        var link = svg.selectAll("path.link")
            .data(links, function (d) {
                return d.target.id;
            });

        // Enter any new links at the parent's previous position.
        link.enter().insert("path", "g")
            .attr("class", "link")
            .attr("d", function (d) {
                var o = {x: source.x0, y: source.y0};
                return diagonal({source: o, target: o});
            });
        // mouseover for edges/paths can be truned off or on
        // .on('mouseover',function(d){
        //     d3.select(this).style("stroke-width",function(d){
        //         return d.target.data.RelAbund*150.0+8;
        //     })
        //     d3.select('#nav_data').html(d3.format(".5f")(d.target.data.RelAbund));
        // })
        // .on('mouseout',function(d){
        //     d3.select(this).style("stroke-width",function(d){
        //         return d.target.data.RelAbund*100.0+4;
        //     })
        //     d3.select('#nav_data').html(d3.format(".5f")(''));
        // });

        // Transition links to their new position.
        link.transition()
            .duration(duration)
            .attr("d", diagonal)
            .style("stroke-width",function(d){
                return d.target.data.RelAbund*100.0+4;
            })
            .style("stroke",function(d){
                if(d.target.data.RelAbund<1e-9){
                    var c = 'rgba(0,0,0,0.0);'
                }else{
                    //console.log(dist_to_color(d.target.data.category,color_range))
                    var c = 'hsla('+dist_to_color(d.target.data.category,color_range)+',100%,50%,0.75)';
                }
                return c;
            });

        // Transition exiting nodes to the parent's new position.
        link.exit().transition()
            .duration(duration)
            .attr("d", function (d) {
                var o = {x: source.x, y: source.y};
                return diagonal({source: o, target: o});
            })
            .remove();

        // Stash the old positions for transition.
        nodes.forEach(function (d) {
            d.x0 = d.x;
            d.y0 = d.y;
        })
    }

    function mouse_over(d){
        var self = d3.select(this);
        self.attr("r", function(d){
            //console.log(d);
            return 60*d.data.RelAbund+8;
        })
            .attr("y", function (d) {
                return d.children || d._children ? -1*(d.data.descendants*30) : 0;
            });
        var display = ' relative abundance for '+d.name;
        d3.select('#nav_data').html(d3.format(".5f")(d.data.RelAbund)+display);
    }

    function mouse_out(d){
        var self = d3.select(this);
        self.attr("r", function (d) {
            return 40 * d.data.RelAbund + 6;
        });
        // if('children' in d || '_children' in d) {
        //     d3.select(this.parentNode).select('text')
        //         .text(function (d) {
        //             return '';
        //         })
        // }
        d3.select('#nav_data').html('');
    }

    function left_click(d){
        if (d.children) {
            d._children = d.children;
            d.children = null;
        }else{
            d.children = d._children;
            d._children = null;
        }
        if(d.data.descendants<1){//leaf click...
            var presence_options = {
                radius: 12,
                fillColor: "hsl(340,100%,50%)",
                color: "hsl(340,100%,50%)",
                opacity: 0.2,
                fillOpacity: 0.6
            };
            var sites = phlyo_to_site_ids(sample_table,d.name);
            var sitesFeature = {type: "FeatureCollection",features:[]};
            for(var i = 0; i < geojsonFeature.features.length; i++){
                if(geojsonFeature.features[i].properties.SID in sites){
                    sitesFeature.features.push(geojsonFeature.features[i]);
                }
            }
            mymap.removeLayer(lastsitesLayer);

            //keep track of the last node
            var node = d3.select(this).node();
            if(node===selected_phis_node) {
                //console.log('already selected...');
                var parent_node = d3.select(this).node().parentNode;
                var parent = d3.select(parent_node).select('circle');
                parent.style('fill',last_phis_color);

                lastsitesLayer = L.geoJSON({type: "FeatureCollection",features:[]},{
                    pointToLayer: function (feature, latlng) {
                        presence_options.radius = 8+30*sites[feature.properties.SID][0];
                        return L.circleMarker(latlng, presence_options);
                    }//,onEachFeature: onEachFeature
                }).addTo(mymap);
                lastsitesLayer.bringToBack();
                selected_phis_node = null;
                update(d);
            }else{ //turn on the layer
                if(selected_phis_node!=null) {
                    var parent_node = d3.select(selected_phis_node).node().parentNode;
                    var parent = d3.select(parent_node).select('circle');
                    parent.style('fill', last_phis_color);
                }

                selected_phis_node = node;
                parent_node = d3.select(this).node().parentNode;
                parent = d3.select(parent_node).select('circle');
                last_phis_color = parent.style('fill');
                parent.style('fill', 'hsla(340,100%,50%,0.6)');

                lastsitesLayer = L.geoJSON(sitesFeature,{
                    pointToLayer: function (feature, latlng) {
                        presence_options.radius = 8+30*sites[feature.properties.SID][0];
                        return L.circleMarker(latlng, presence_options);
                    }//,onEachFeature: onEachFeature
                }).addTo(mymap);
                lastsitesLayer.bringToBack();
            }
            //console.log(selected_phis_node);
        }else{
            update(d);
        }
    }
    window.addEventListener("resize",update);
}

//build a auto DOM
$(function() {
    //upper data bar for display of SID, relative abundance and species name
    $('<div id="nav_bar"></div>').appendTo('body');
    $('<span id="nav_site">phis_viz</span>').appendTo('#nav_bar');
    $('<span id="nav_data"></span>').appendTo('#nav_bar');
    $('<span id="map_toggle"><button id="map_toggle_button">map\></button></span>').appendTo('#nav_bar');

    //lower temporal bar for drilling down into a SID data set once it has already been selected
    $('<div id="temporal_bar"></div>').appendTo('body')

    //left and right graph panel containers
    $('<div id="phis_viz"></div>').appendTo('body');
    $('<span id="phis_side"></span>').appendTo('#phis_viz');
    $('<span id="map_side"></span>').appendTo('#phis_viz');
    phlyo_tree_map_graph('data/phylo_tree.csv','data/sample_data.csv','phis_side',
                         'data/site_geo.json','map_side',true);

});