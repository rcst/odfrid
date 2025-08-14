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

SELECT trip_id, ton.route_id, arrival_time, departure_time, sequence, 
CASE 
  WHEN direction = 'F' THEN false
  ELSE true 
END AS inbound,
ST_SetSRID(ST_MakePoint(stop_lon, stop_lat), 4326) AS geom
FROM survey_on_bus_stop AS son
JOIN survey_on_bus_trip AS ton USING (trip_id)
WHERE son.route_id = 2

