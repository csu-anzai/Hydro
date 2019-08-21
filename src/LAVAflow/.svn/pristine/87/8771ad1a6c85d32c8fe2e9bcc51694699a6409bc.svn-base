function weight = WeightFinder(P,Q,R, PR,PQ,RQ)

	Qa = acos(dot(PQ,RQ)/(norm(PQ,2)*norm(RQ,2)));
	Ra = acos(dot(PR,RQ)/(norm(PR,2)*norm(RQ,2)));
	Pa = acos(dot(PR,PQ)/(norm(PR,2)*norm(PQ,2)));


	if Pa <= pi/2 & Qa <= pi/2 & Ra <= pi/2

		weight = 1/8*(norm(PR,2)^2*cot(Qa)+ norm(PQ,2)^2*cot(Ra));

	else if Pa <= pi/2
		weight = 1/8 * norm(cross(PQ,PR),2);
	else 
 		weight = 1/4 * norm(cross(PQ,PR),2);

	end
end
