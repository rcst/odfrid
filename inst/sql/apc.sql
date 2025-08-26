-- SELECT routes.route_number, stops.name, son.*
-- FROM survey_on_bus_stop AS son
-- JOIN routes USING (route_id)
-- JOIN stops USING (stop_id)
-- WHERE route_number = '2'
-- ORDER BY arrival_time;

-- SELECT trip_id, COUNT(*) AS nstops
-- FROM survey_on_bus_stop AS son
-- JOIN routes USING (route_id)
-- JOIN stops USING (stop_id)
-- WHERE route_number = '2'
-- GROUP BY 1

--WITH cte AS (
--	SELECT trip_id, COUNT(*) AS nstops
--	FROM survey_on_bus_stop AS son
--	JOIN routes USING (route_id)
--	JOIN stops USING (stop_id)
--	WHERE route_number = '2'
--	GROUP BY 1
--)
--SELECT nstops, COUNT(*) 
--FROM cte
--GROUP BY 1
--ORDER BY nstops;

-- SELECT *
-- FROM route_stops
-- WHERE route_id = 2
-- ORDER BY seq_id

-- example 13 stops

-- SELECT stops.name, son.* 
-- FROM survey_on_bus_stop AS son
-- JOIN stops USING (stop_id)
-- WHERE trip_id = 38671153
-- ORDER BY arrival_time

WITH cte AS (
	WITH cte AS (
		SELECT trip_id, ton.route_id, arrival_time, departure_time, sequence,	  
		CASE	  
			WHEN direction = 'F' THEN false	 
			ELSE true	  
	END AS inbound,	 
	ST_SetSRID(ST_MakePoint(stop_lon, stop_lat), 4326) AS geom	 
	FROM survey_on_bus_stop AS son	 
	JOIN survey_on_bus_trip AS ton USING (trip_id) 
) 
SELECT * 
FROM cte 
JOIN LATERAL (	 
	SELECT stop_id, rs.seq_id, s.geom AS mapped_geom, ST_DISTANCE(s.geom, ST_TRANSFORM(cte.geom, 32634)) AS dist_mapped_geom
	FROM route_stops AS rs	 
	JOIN stops AS s USING (stop_id)	 
	WHERE rs.route_id = cte.route_id	 
	AND rs.inbound = cte.inbound	 
	ORDER BY ST_Transform(s.geom, 4326) <-> cte.geom	 
	LIMIT 1 ) 
AS rs ON true
)
SELECT trip_id, route_id, arrival_time, departure_time, sequence, seq_id, dist_mapped_geom AS d, inbound
FROM cte
WHERE route_id = 2
AND inbound = false;
