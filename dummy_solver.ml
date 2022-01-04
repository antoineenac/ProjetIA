open Modules
   
let conflict dmin =
  let dmin_diag= dmin *. sqrt 2. in
  let epsilon= 1e-5 in
  let dbox= dmin_diag+.epsilon in
  fun pj pk ->
    let in_box =
      (abs_float (pj.Geo.P2D.x-.pk.Geo.P2D.x))<= dbox
      && (abs_float (pj.Geo.P2D.y-.pk.Geo.P2D.y))<= dbox in
    in_box && Geo.P2D.distance pj pk<=dmin

(*------------------------------*)

let avoid dmin others f =
  (* Move [f] away from neighbors *)
  let c= {Geo.V2D.x=0.;y=0.} in
  let pos= F.position f in
  List.fold_left
    (fun c other ->
      let p= F.position other in
      if conflict dmin p pos then
        let open Geo.V2D in
        c + Geo.V2D.make p pos
      else c)
    c others                     

let to_destination pln f =
  Geo.V2D.make (F.position f) (P.arr pln)

(*------------------------------*)

let closest_admissible_direction (lb,ub) dir =
  if Arc.inside (lb,ub) dir then dir
  else
    let alpha= Arc.closest_angular_distance dir lb
    and beta = Arc.closest_angular_distance dir ub in
    if alpha <= beta then lb else ub

(*------------------------------*)

let c1= 1.;; (* Avoid conflicts *)
let c2= 1e-2;; (* Go to destination *)

let compute_velocities max_turn_angle dmin plns flying =
  let rec loop prev l acc =
    match l with
      [] -> acc
    | f::ls ->
       let current_dir= F.vdir f in
       let lb= current_dir -. max_turn_angle
       and ub= current_dir +. max_turn_angle in
       let others= List.rev_append prev ls in
       let v_avoid= avoid dmin others f in
       try
         let pln= P.find (F.id f) plns in
         let v_to_dest= to_destination pln f in
         (*let new_v= Geo.V2D.( *~ ) c2 v_to_dest in *)
         let new_v=
           Geo.V2D.(+) (Geo.V2D.( *~ ) c1 v_avoid)
             (Geo.V2D.( *~ ) c2 v_to_dest) in
         let new_vdir= Geo.V2D.direction new_v in
         let new_vdir= closest_admissible_direction (lb,ub) new_vdir in
         let vnorm= F.vnorm f in         
         let new_acc= Util.IntMap.add (F.id f) (new_vdir, vnorm) acc in
         loop (f::prev) ls new_acc
       with Not_found -> failwith "compute_velocities" in
  loop [] flying Util.IntMap.empty

(*------------------------------------------------------*)
 
let deep_copy_pop = fun pop ->
	let taille_pop = Array.length pop in
	let taille_individu = Array.length pop.(0) in
	let copy = Array.make_matrix taille_pop taille_individu pop.(0).(0) in
	for k = 0 to taille_pop-1 do
		for i = 0 to taille_individu-1 do
			copy.(k).(i)<- pop.(k).(i)
		done;
	done;
	copy

open Random;;

Random.self_init;;

let signe = fun x -> if x < 0. then - 1. else 1.;;
    

let one_step_evol_diff_float = fun varmax pop fobj cr f dirobj ->
    let nb_individu = Array.length pop in
    let taille_individu = Array.length pop.(0) in
    let chosen = Array.make 3 0 in
    for k=0 to (nb_individu-1) do
		let trial = Array.make taille_individu 0.0 in
		chosen.(0)<-Random.int (nb_individu-1);
      chosen.(1)<-Random.int (nb_individu-1);
      chosen.(2)<-Random.int (nb_individu-1);
      
		while chosen.(0)=k do 
			chosen.(0)<- Random.int (nb_individu-1);	
	   done;
	   
	   (*Printf.printf "\n \n k = %d \n" k;
	   Printf.printf "C0 = %d \n" (chosen.(0));*)
	   
      while chosen.(1)=chosen.(0) || chosen.(1)=k do 
			chosen.(1)<-Random.int (nb_individu-1);
	   done;
	   
	   (*Printf.printf "C1 = %d \n" (chosen.(1));*)
	   
	   while chosen.(2) = chosen.(0) || chosen.(2) = chosen.(1) || chosen.(2) = k do
			chosen.(2) <- Random.int (nb_individu-1)
		done;
		
	   (*Printf.printf "C2 = %d \n" (chosen.(2));*)
	   
      let a = chosen.(0) in 
	   let b = chosen.(1) in 
	   let c = chosen.(2) in
      for i=0 to (taille_individu-1) do
			if (Random.int 100)<cr then
				let change = (pop.(a).(i) +. f*.(pop.(b).(i)-.pop.(c).(i)) -. pop.(k).(i)) in
				if abs_float(change) > varmax 
					then 
						trial.(i) <- pop.(k).(i) +. varmax*.signe(change)
					else  
						trial.(i) <- pop.(a).(i) +. f*.(pop.(b).(i)-.pop.(c).(i))
         else
            trial.(i) <- pop.(k).(i)
    done;
    if (fobj trial dirobj)<(fobj pop.(k) dirobj) then pop.(k) <- trial
    done;
    pop;;
    
let evol_diff_float = fun varmax pop gen_max fobj cr f dirobj ->
    let aux_pop = ref (deep_copy_pop pop) in
    for k=0 to gen_max do
        aux_pop := one_step_evol_diff_float varmax (!aux_pop) fobj cr f dirobj;
    done;
    !aux_pop;;
    
let fobj trial dirobj = 
	let n = Array.length dirobj in 
	let ecart = ref 0. in 
	for i=0 to (n-1) do 
		ecart := !ecart +. (abs_float (trial.(i) -. dirobj.(i)))**2.
	done; 
	!ecart
	
let gather_dir= fun flying ->
	let n = List.length flying in 
	
	let dir_vect = Array.make n 0. in 
	let rec aux flying i= 
		match flying with
		|[]-> dir_vect
		|f::tl-> dir_vect.(i)<-F.vdir f;
			 aux tl (i+1) in 
	aux flying 0
	
let gather_pos= fun flying ->
	let n = List.length flying in 
	let pos_vect = Array.make n {Geo.P2D.x=0.;y=0.} in
	let rec aux flying i= 
		match flying with
		|[]-> pos_vect
		|f::tl-> pos_vect.(i)<-F.position f;
			 aux tl (i+1) in 
	aux flying 0

let gather_norm= fun flying ->
	let n = List.length flying in 
	let norm_vect = Array.make n 0. in
	let rec aux flying i= 
		match flying with
		|[]-> norm_vect
		|f::tl-> norm_vect.(i)<-F.vnorm f;
			 aux tl (i+1) in 
	aux flying 0
	
let gather_dirobj= fun flying plns ->
	let n = List.length flying in 
	let dirobj_vect = Array.make n 0. in
	let rec aux flying i= 
		match flying with
		|[]-> dirobj_vect
		|f::tl-> let pln= P.find (F.id f) plns in
         let v_to_dest= to_destination pln f in dirobj_vect.(i)<-(Geo.V2D.direction v_to_dest);
			 aux tl (i+1) in 
	aux flying 0
	
let make_pop_test_float = fun taille_pop taille_individu ->
	let pop = Array.make_matrix taille_pop taille_individu 0.0 in
	for k=0 to taille_pop-1 do
		for i=0 to taille_individu-1 do
			let aux = Random.float 6.28 in 
			pop.(k).(i) <- aux;
		done;
	done;
	pop;;


let compute_velocities_2 max_turn_angle plns flying =
  let n = List.length flying in 
  let n_pop = 20 in
  let n_gen = 500 in
  let cr = 66 in 
  let f = 0.6 in  
	let dir_vect = gather_dir flying in 
	let norm_vect = gather_norm flying in 
	let dirobj_vect = gather_dirobj flying plns in 
	if n>0 then 
		let pop = make_pop_test_float n_pop n in 
  	let new_v = evol_diff_float max_turn_angle pop n_gen fobj cr f dirobj_vect in 
  	let rec loop l acc i  =
  	match l with 
  	| [] -> acc
  	| f:: tl -> 
  		try
  			let new_dir = new_v.(0).(i) in 
  			let new_norm = norm_vect.(i) in 
  			let new_acc= Util.IntMap.add (F.id f) (new_dir, new_norm) acc in
  			loop tl new_acc (i+1)
  		with Not_found -> failwith "compute_velocities" in
		loop flying Util.IntMap.empty 0
	else 
		let rec loop l acc i  =
  	match l with 
  	| [] -> acc
  	| f:: tl -> 
  		try
  			let new_dir = dir_vect.(i) in 
  			let new_norm = norm_vect.(i) in 
  			let new_acc= Util.IntMap.add (F.id f) (new_dir, new_norm) acc in
  			loop tl new_acc (i+1)
  		with Not_found -> failwith "compute_velocities" in
		loop flying Util.IntMap.empty 0
			
		
		
		
