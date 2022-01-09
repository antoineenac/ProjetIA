(***********************************************************************)
(*                                                                     *)
(*             Smallsimu : a small air traffic simulation              *)
(*                                                                     *)
(*         David Gianazza, Ecole Nationale de l'Aviation Civile        *)
(*                                                                     *)
(*  Copyright 2020 Ecole Nationale de l'Aviation Civile.               *)
(*  All rights reserved.  This file is distributed under the terms of  *)
(*  the GNU Library General Public License.                            *)
(*                                                                     *)
(***********************************************************************)

include Modules

let debug= ref false
   
type user_param= {
    initial_plns: P.db;
    twindow: (Sim_types.second * Sim_types.second) option;
    (* Time interval used in the performance metrics calculation *)
    max_turning_rate: Sim_types.radian_per_second;
    selected_flights: Util.IntSet.t option;
    debug_poisson: (Sim_types.meter * Pb.debug_map) option;
    plotfile: string }

(*
type param=
  { ti: Sim_types.second;
    tf: Sim_types.second;
    timestep: Sim_types.second;
    lateral_sep: Sim_types.meter;
    sep_at_dep: Sim_types.meter;
    tma_radius: Sim_types.meter;
    anticipation: Sim_types.second; (* parameter tau, in ORCA *)
  }
 *)
               
open Ustate
                      
(*------------------------------*)
   
let solve uparam param state flying =
(*
  let velocities=
    F.fold
      (fun id f acc ->
        Util.IntMap.add id (F.vdir f, F.vnorm f) (*F.velocity f*) acc)
  flying Util.IntMap.empty in
 *)
  let dmin= 10. *. Util.nm2meter in
  let lf= List.map snd (F.list flying) in
  Printf.printf "Flying flights:\n";
  List.iter (fun f -> Printf.printf " %d" (F.id f)) lf;
  Printf.printf "\n";flush stdout;
  let max_turn_angle= param.Sim.timestep*.uparam.max_turning_rate in
  let velocities=
    Dummy_solver.compute_velocities max_turn_angle dmin state.plns lf 
    in
  velocities, ()

(*------------------------------*)
   
let do_at_init uparam _param _plns = ()

let do_at_iter uparam _param state _history =
  let nf= F.cardinal state.flights in
  let nc= List.length state.sep_losses in
  Printf.printf "%d:  t= %.2f  nf= %d  nc=%d  " state.iter state.time nf nc;
  flush stdout
(*Geo.debug_flag:= (k= 1071);*)
(*    if k= 1070 then Geo.debug_flag:= true; *)  
    
let do_after_solving uparam _param state _history =
  let nc= List.length state.sep_losses in
  Printf.printf "remaining= %d\n" nc;flush stdout
  
let do_at_exit uparam param state history =
  Printf.printf "\n";
  if state.time>=param.Sim.tf then (
    Printf.printf "Exit at t>=%.2f\n" param.tf;flush stdout);
  Printf.printf "ti= %.2f   tf= %.2f\n" param.ti param.tf;
  let (t1,t2) =
    match uparam.twindow with
    | None -> (neg_infinity,infinity)
    | Some (t1,t2) ->
       Printf.printf "Metrics measured between %.2f and %.2f\n" t1 t2;
       flush stdout;
       (t1,t2) in
  H.save_plots None history uparam.plotfile;
  if not !Options.display then (
    let cmd=
      Printf.sprintf "gnuplot -p -e \"plot \'%s\' w l\"" uparam.plotfile in
    let exit_code= Sys.command cmd in
    Printf.printf "Finished %d\n" exit_code);
  let res=
    H.remaining_separation_losses ~start:t1 ~stop:t2 uparam.selected_flights
      param.lateral_sep history in
  Printf.printf "Number of pointwise separation losses: %d\n"
    res.H.nb_pointwise_conflicts;
  Printf.printf "Number of sep. losses between flights: %d\n"
    res.H.nb_flight_conflicts;
  Printf.printf "Number of sep. losses between trajectory segments: %d\n"
    res.H.nb_traj_conflicts;
  Printf.printf "Total duration of separation losses (in seconds): %.2f\n"
    res.H.total_duration;
  flush stdout;
  let plns= H.Ustate.plns state in
  let trajectories= H.trajectory_map uparam.selected_flights history in
  let (terminated,active)= H.terminated_and_active trajectories in
  let lengthening =
    (H.total_traj_lengthening ~start:t1 ~stop:t2 terminated)/.Util.nm2meter in
  let lengthening2 =
    (H.total_traj_lengthening2 ~start:t1 ~stop:t2 uparam.initial_plns
       trajectories)
    /.Util.nm2meter in
  let total_delay= H.total_delay ~start:t1 ~stop:t2 trajectories in
  let not_terminated= Util.IntMap.cardinal active in                        
  let total_remaining_dist=
    (H.total_remaining_distance ~start:t1 ~stop:t2 plns active)
    /.Util.nm2meter in
  Printf.printf "Total number of flights: %d\n"
    (Util.IntMap.cardinal trajectories);
  Printf.printf "Total trajectory lengthening (in NM) for terminated flights: %.2f\n" lengthening;
  Printf.printf "Total trajectory lengthening (in NM) for all flights (counting the shortest distance from last position to arrival for flights still active) : %.2f\n" lengthening2;    
  Printf.printf "Total delays (in seconds) for all flights: %.2f\n"
    total_delay;
  Printf.printf "Number of flights still active: %d\n" not_terminated;
  Printf.printf "Total remaining distance (in NM) for flights still active: %.2f\n" total_remaining_dist;
  Printf.printf "FINAL SCORE: (%f,%f,%f)\n" res.H.total_duration lengthening2 total_delay;
  flush stdout;
  H.save_flight_counts history "flight_counts"
    
(*------------------------------*)
  
type scaling_fun= (float*float -> int*int)
         
let departing state history=
  match history with
  | prev::_->
     F.list
       (F.filter
          (fun id f ->
            (F.is_active f || F.is_delayed f) &&
              (try let _f2= F.find id prev.Ustate.flights in false
               with Not_found -> true ) )
          state.Ustate.flights)
  | [] -> F.list state.Ustate.flights

let debug_tag= "debug"

let draw_exclusion_zones cv scale exclusion_zones =
  let tags= [Tkview.scalable_tag;Tkview.flight_tag;debug_tag] in
  let scale3 {Util.SmallGeo.x;y} = scale (x,y) in
  let color= `Color "purple" in
  List.iter
    (fun ({Util.SmallGeo.c;r} as circle) ->
      let (p1,p2)= Util.SmallGeo.circle_bbox circle in
      let (x1,y1)= scale3 p1 and (x2,y2)= scale3 p2 in
      ignore (Canvas.create_oval ~x1:x1 ~y1:y1 ~x2:x2 ~y2:y2
                ~outline:color ~width:1
                ~tags:tags cv) )
    exclusion_zones
  
let debug_print= ref false
               
let draw_domain cv scale color radius domain =
  let tags= [Tkview.scalable_tag;Tkview.flight_tag;debug_tag] in
  let scale3 {Util.SmallGeo.x;y} = scale (x,y) in
  let big_circle= {Util.SmallGeo.c= {Util.SmallGeo.x=0.;y=0.};r=radius} in
  let (p1,p2)= Util.SmallGeo.circle_bbox big_circle in
  let (x1,y1)= scale3 p1 and (x2,y2)= scale3 p2 in
  List.iter
    (fun (lb,ub) ->
      let lb_deg= Util.rad2deg (Arc.normalize_angle 0. lb) in
      let extent= Arc.trigo_angle_diff lb ub in
      let extent_deg= Util.rad2deg extent in
      if !debug_print then (
        Printf.printf "lb= %.2f  extent= %.2f\n" lb_deg extent_deg;
        flush stdout);
      ignore (Canvas.create_arc ~x1:x1 ~y1:y1 ~x2:x2 ~y2:y2
                ~start:lb_deg ~extent:extent_deg
                ~outline:color
                ~width:2
                (*~fill:arc_color*)
                (*~stipple: (`Predefined "gray50")*)
                ~tags:tags cv) )
    domain
  
let draw_ids cv scale state =
  let dx= 10 and dy=10 in
  let scale2 p = scale (G.x p, G.y p) in
  F.iter
    (fun id f ->
      if F.is_active f || F.is_delayed f then (
        let tags= [Tkview.scalable_tag;Tkview.flight_tag;debug_tag(*;flight_id_tag f*)] in
        let p= F.position f in
        let x,y= scale2 p in
        let txt= Printf.sprintf "%d" id in
        ignore (Canvas.create_text ~x:(x+dx) ~y:(y+dy)
                  ~text:txt ~justify:`Left
                  ~width:50
                  ~tags:tags cv) ) )
    state.Ustate.flights
  
let debug_drawing uparam param view scale state history =
  draw_ids view.Tkview.cv scale state;;

let do_my_drawing user_param param view scale state history =
  let cv= view.Tkview.cv in
  if !debug then
    debug_drawing user_param param view scale state history;
  let domain_color= `Color "green"
  and forbidden_color= `Color "red" in
  match user_param.debug_poisson with
  | None -> ()
  | Some (radius,map) ->
     let l= departing state history in
     match l with
       [] ->
        Printf.printf "DELETING DEBUG DATA\n";flush stdout;
        Canvas.delete cv [`Tag debug_tag]
     | _->
        List.iter
          (fun (id,pln) ->
            try
              let {Pb.domain;forbidden;exclusion_zones}=
                Util.IntMap.find id map in
              Printf.printf "DOMAIN:\n";
              Arc.print_arcs domain;flush stdout;
              Printf.printf "\nFORBIDDEN:\n";
              Arc.print_arcs forbidden;
              Printf.printf "\n";flush stdout;
              draw_exclusion_zones cv scale exclusion_zones;
              debug_print:=true;
              draw_domain cv scale forbidden_color radius forbidden;
              debug_print:=false;
              draw_domain cv scale domain_color radius domain
            with Not_found ->
              failwith (Printf.sprintf "do_my_drawing : %d not found\n" id) )
          l
(*          raise Util.Debug *)
       
let my_buttons uparam param top view =
  let c_frame= Frame.create top in
  let relief= if !debug then `Sunken else `Flat in
  let btn_debug = Button.create ~text:"Debug" ~relief:relief c_frame in      
  Tk.pack ~expand:true ~fill:`Both [c_frame];
  Tk.pack ~expand:true ~fill:`X ~side:`Left [btn_debug];
  [|btn_debug|]
  
let configure_my_buttons uparam param _top view scale state history =
  let cmd_debug = fun ()->
    Printf.printf "Debug\n";flush stdout;
    debug:= not !debug;
    let relief= if !debug then `Sunken else `Flat in
    Button.configure ~relief:relief view.Tkview.my_buttons.(0);
    Canvas.delete view.Tkview.cv [`Tag debug_tag];
    do_my_drawing uparam param view scale state history in
  Button.configure ~command:cmd_debug view.Tkview.my_buttons.(0)

