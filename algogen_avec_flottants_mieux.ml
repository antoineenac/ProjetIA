open Random;;

Random.self_init;;


let deep_copy_pop = fun pop ->
	let taille_pop = Array.length pop in
	let taille_individu = Array.length pop.(0) in
	let copy = Array.make_matrix taille_pop taille_individu pop.(0).(0) in
	for k = 0 to taille_pop-1 do
		for i = 0 to taille_individu-1 do
			copy.(k).(i)<- pop.(k).(i)
		done;
	done;
	copy;;


	
let close_zero_float = fun tab ->
	let m = ref 0.0 in
	for i=0 to Array.length tab-1 do
		if tab.(i)< 0.0 then m:= !m -. tab.(i) else m:= !m +. tab.(i)
	done;
	!m ;;


let signe = fun x -> if x < 0. then - 1. else 1.;;
    

let one_step_evol_diff_float = fun varmax vect_precedent pop fobj cr f ->
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
				let change = (pop.(a).(i) +. f*.(pop.(b).(i)-.pop.(c).(i)) -. vect_precedent.(i)) in
				if abs_float(change) > varmax 
					then 
						trial.(i) <- vect_precedent.(i) +. varmax*.signe(change)
					else  
						trial.(i) <- pop.(a).(i) +. f*.(pop.(b).(i)-.pop.(c).(i))
         else
            trial.(i) <- pop.(k).(i)
    done;
    if (fobj trial)<(fobj pop.(k)) then pop.(k) <- trial
    done;
    pop;;

    
let evol_diff_float = fun varmax vect_precedent pop gen_max fobj cr f ->
    let aux_pop = ref (deep_copy_pop pop) in
    for k=0 to gen_max do
        aux_pop := one_step_evol_diff_float varmax vect_precedent(!aux_pop) fobj cr f;
    done;
    !aux_pop;;

let make_pop_test_float = fun taille_pop taille_individu ->
	let pop = Array.make_matrix taille_pop taille_individu 0.0 in
	for k=0 to taille_pop-1 do
		for i=0 to taille_individu-1 do
			let aux = Random.float 360. in 
			pop.(k).(i) <- aux;
		done;
	done;
	pop;;


let varmax = 3.;;

let pop_test_float = make_pop_test_float 10 10;;

let vect_precedent = [|250.;10.;64.;80.;90.;1.;190.;300.;64.;20.|];;

let test = evol_diff_float varmax vect_precedent pop_test_float 1000 close_zero_float 53 0.6;;




let best_individu = fun pop fobj ->
	let nb_individu = Array.length pop in 
	let min = ref (fobj pop.(0)) in
	let indice_min = ref 0 in 
	for k=1 to nb_individu-1 do 
		let aux = (fobj pop.(k)) in
		if aux <(!min) then begin min := aux; indice_min := k end
	done;
	!indice_min,!min;;
	
let best_test = best_individu test close_zero_float;;

let optim_cr_f = fun nb_gen taille_pop taille_individu fobj ->
	let pop = ref (make_pop_test_float taille_pop taille_individu) in
	let score = ref infinity in
	let cr = ref 0 in
	let bestcr = ref 0 in 
	let f = ref 0.0 in
	let bestf = ref 0.0 in 
	for pas_cr = 0 to 100 do 
		cr := !cr + 1 ;
		f := 0.0;
		for pas_f = 0 to 20 do
			f := !f +. 0.1;
			let aux = evol_diff_float !pop nb_gen fobj !cr !f in 
			let couple = best_individu aux fobj in
			let meilleur_score = snd couple in
			if meilleur_score < !score then begin
				score := meilleur_score;
				bestcr := !cr;
				bestf := !f;
			end
		done;
		Printf.printf "%d \n" !cr;
	done;
	!bestcr,!bestf,!score;;

let test_optim = optim_cr_f 1000 10 10 close_zero_float;;			

let occurence_optim = fun nb_test nb_gen taille_pop taille_individu fobj ->
	let occ = Array.make nb_test (0,0.0,0.0) in
	for k=0 to nb_test-1 do
		occ.(k)<- (optim_cr_f nb_gen taille_pop taille_individu fobj)
	done;
	occ;;

(*let occ = occurence_optim 20 1000 10 5 close_zero_float;;*)
(*		[|(66, 0.6, 1.1900320215263738e-047); (77, 0.6, 8.2427788905984683e-047);
    (63, 0.6, 6.6006649311369449e-046); (68, 0.6, 1.0051197117003727e-045);
    (68, 0.6, 6.0786206306978488e-048); (64, 0.6, 2.4581558653121063e-045);
    (68, 0.6, 4.8677760380082514e-048); (78, 0.6, 2.3908344848862848e-046);
    (56, 0.6, 2.6491611588254992e-044); (67, 0.6, 8.9080721802556135e-046);
    (50, 0.5, 2.2616578254220583e-047); (61, 0.6, 2.3402160246592014e-045);
    (48, 0.5, 3.0325959029345727e-045); (65, 0.6, 5.0328410339730423e-046);
    (70, 0.6, 4.4723107653625039e-045); (60, 0.6, 2.2701851822439157e-047);
    (65, 0.6, 7.4883295051860618e-046); (71, 0.6, 8.923481595077688e-046);
    (66, 0.6, 3.0965192944897027e-046); (66, 0.6, 3.3216596686122976e-046)|] *)
		

