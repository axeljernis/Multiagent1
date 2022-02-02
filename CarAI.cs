using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.Linq;


namespace UnityStandardAssets.Vehicles.Car
{
    [RequireComponent(typeof(CarController))]
    public class CarAI : MonoBehaviour
    {
        private CarController m_Car; // the car controller we want to use
        private float[] carsize = new float[] { 3, 5 };

        public GameObject terrain_manager_game_object;
        TerrainManager terrain_manager;

        List<Vector3> my_path;
        public List<float> pathAngles;

        public List<Vector3> path = new List<Vector3>();


        private void Start()
        {
            // get the car controller
            m_Car = GetComponent<CarController>();
            terrain_manager = terrain_manager_game_object.GetComponent<TerrainManager>();

            // Plan your path here
            // Replace the code below that makes a random path
            // ...

            Vector3 start_pos = terrain_manager.myInfo.start_pos;
            Vector3 goal_pos = terrain_manager.myInfo.goal_pos;
            TerrainInfo tinfo = terrain_manager.myInfo;
            float[,] newMap = traversableMap(tinfo.traversability);
            newMap = InflateObstacles(newMap, carsize, tinfo);

            my_path = this.SearchPath(start_pos, goal_pos, newMap, tinfo);

            //my_path[0] = start_pos;
            //my_path[1] = (m_Car.transform.forward + start_pos);


            pathAngles = new List<float>();
            pathAngles.Add(0f);

            //Create angles at all intersections
            for (int i = 0; i < my_path.Count-2; i++)
            {
                my_path[i] = CheckIfWaypointcollidesAndAdjust(my_path[i], newMap, tinfo, carsize);
                float angle = PathSegment.getAngleForMiddle(my_path[i], my_path[i+2]);
                pathAngles.Add(angle);
                //PathSegment.createVertexLine(angle - (float)(Math.PI / 2), my_path[1]);
            }


            /*
            float angle = PathSegment.getAngleForMiddle(my_path[1], my_path[2]);
            float cutoffAngle = PathSegment.getAngleForMiddle(my_path[0], my_path[1]);

            PathSegment.createPathSeg(0, (float)Math.PI / 2, my_path[0], my_path[1], cutoffAngle);


            float oldAngle = angle;
            angle = PathSegment.getAngleForMiddle(my_path[3], my_path[2]);
            //cutoffAngle = PathSegment.getAngleForMiddle(my_path[1], my_path[2]);

            PathSegment.createPathSeg((float)Math.PI / 2, angle, my_path[1], my_path[2], (float)Math.PI / 2);
            */



            // Plot your path to see if it makes sense
            // Note that path can only be seen in "Scene" window, not "Game" window
            Vector3 old_wp = start_pos;
            foreach (var wp in my_path)
            {
                Debug.DrawLine(old_wp, wp, Color.red, 100f);
                old_wp = wp;
            }

            
        }

        private List<Vector3> SearchPath(Vector3 start_pos, Vector3 goal_pos, float[,] newMap, TerrainInfo tinfo)
        {
            this.path.Clear();

            //Directions to neighbouring cells
            List<Vector2> neighbours = new List<Vector2>() {
                new Vector2(0,3),
                new Vector2(3,0),
                new Vector2(0,-3),
                new Vector2(-3,0),
                new Vector2(2,2),
                new Vector2(-2,-2),
                new Vector2(2,-2),
                new Vector2(-2,2),

                new Vector2(1,2),
                new Vector2(2,1),

                new Vector2(-1,-2),
                new Vector2(-2,-1),

                new Vector2(1,-2),
                new Vector2(2,-1),

                new Vector2(-1,2),
                new Vector2(-2,1),


                new Vector2(1,3),
                new Vector2(3,1),

                new Vector2(-1,-3),
                new Vector2(-3,-1),

                new Vector2(1,-3),
                new Vector2(3,-1),

                new Vector2(-1,3),
                new Vector2(-3,1),

                new Vector2(2,3),
                new Vector2(3,2),

                new Vector2(-2,-3),
                new Vector2(-3,-2),

                new Vector2(2,-3),
                new Vector2(3,-2),

                new Vector2(-2,3),
                new Vector2(-3,2),


            };

            List<PathPoint> open = new List<PathPoint>();
            List<PathPoint> closed = new List<PathPoint>();
            List<PathPoint> temp_path = new List<PathPoint>();


            //Get the starting position for the A* algorithm
            int i = large_get_i_index(start_pos.x, newMap, tinfo);
            int j = large_get_j_index(start_pos.z, newMap, tinfo);

            PathPoint startingPoint = new PathPoint();
            startingPoint.SetParent(startingPoint);
            startingPoint.i = i;
            startingPoint.j = j;
            startingPoint.h = 0f;
            startingPoint.g = 0f;
            startingPoint.f = 0f;
            startingPoint.direction = new Vector2(0f, 1f);

            open.Add(startingPoint);

            //Get the goal position for the A* algorithm
            int i_goal = large_get_i_index(goal_pos.x, newMap, tinfo);
            int j_goal = large_get_j_index(goal_pos.z, newMap, tinfo);

            PathPoint lastPos = startingPoint;
            Debug.Log(m_Car.transform.forward.x + " " + m_Car.transform.forward.y);
            int counter = 0;
            //while last_pos not goal position
            while ((lastPos.i != i_goal || lastPos.j != j_goal) && counter <= 2000000)
            {
                counter += 1;
                foreach (Vector2 direction in neighbours)
                {
                    int new_i = (int)(lastPos.i + direction.x);
                    int new_j = (int)(lastPos.j + direction.y);

                    //Check if the neighbour is on the map
                    if (newMap[new_i, new_j] == null) continue;
                    //Check if the neighbour is traversable
                    if (newMap[new_i, new_j] > 0.5)
                    {
                        continue;
                    }

                    //Create A PathPoint object for this location
                    PathPoint pathP = new PathPoint();
                    pathP.SetParent(lastPos);
                    pathP.i = new_i;
                    pathP.j = new_j;

                    //Check if the neighbour is already in the closed list
                    bool is_closed = false;
                    foreach (PathPoint p in closed)
                    {
                        if (p.i == pathP.i && p.j == pathP.j)
                        {
                            is_closed = true;
                        }
                    }
                    if (is_closed) continue;



                    List<PathPoint> currentPath = new List<PathPoint>();
                    PathPoint currentNode = pathP;
                    int pathLength = 0;
                    //while (currentNode != startingPoint)
                    //{
                    currentPath.Add(currentNode);
                    currentNode = currentNode.parentNode;
                    pathLength += 1;
                    //}
                    //Debug.Log(pathLength);
                    Vector2 lastdirection;
                    if (temp_path.Count == 10000000000000)
                    {
                        lastdirection = m_Car.transform.forward;
                    }
                    else
                    {
                        lastdirection = lastPos.direction;
                    }

                    float anglediff = Vector2.Angle(lastdirection, direction);
                    Debug.Log(anglediff + "angleddddddddif" + lastdirection.x + " " + lastdirection.y + " " + direction.x);
                    //calculate the distance travelled (distance to lastpos + lastPos distance to start)
                    float g = Vector2.Distance(new Vector2(lastPos.i, lastPos.j), new Vector2(new_i, new_j)) + lastPos.g;
                    //calculate the distance to the goal
                    float h = Vector2.Distance(new Vector2(new_i, new_j), new Vector2(large_get_i_index(goal_pos.x,newMap, tinfo), large_get_j_index(goal_pos.z, newMap, tinfo))) + (anglediff*anglediff);//add expression for angle diff
                    float f = g + h;

                    //Add the distance and fitness values
                    pathP.h = h;
                    pathP.g = g;
                    pathP.f = f;
                    pathP.direction = direction;
                    //add it to the open list

                    bool is_open = false;
                    foreach (PathPoint p in open)
                    {
                        if ((p.i == pathP.i && p.j == pathP.j))
                        {
                            if (p.f > pathP.f)
                            {
                                p.parentNode = pathP.parentNode;
                                p.f = pathP.f;
                                p.direction = pathP.direction;
                            }
                            is_open = true;
                            //Debug.Log("already in open");
                        }
                    }
                    if (is_open) continue;
                    open.Add(pathP);

                }
                closed.Add(lastPos);
                open.Remove(lastPos);

                //order the open list by its f value
                open = open.OrderBy(p => p.f).ToList<PathPoint>();
                lastPos = open.ElementAt(0);
                temp_path.Add(lastPos);
            }
            if (lastPos.i == i_goal && lastPos.j == j_goal)
            {
                Debug.Log("Path found");
            }
            temp_path = this.RetracePath(startingPoint, lastPos, temp_path);
            foreach (PathPoint p in temp_path)
            {
                path.Add(new Vector3(large_get_i_pos(p.i, newMap, tinfo), 0, large_get_j_pos(p.j, newMap, tinfo)));
            }
            return path;
        }

        List<PathPoint> RetracePath(PathPoint startNode, PathPoint targetNode, List<PathPoint> path)
        {
            List<PathPoint> retracedPath = new List<PathPoint>();
            PathPoint currentNode = targetNode;


            while (currentNode != startNode)
            {
                retracedPath.Add(currentNode);
                currentNode = currentNode.parentNode;
            }
            //retracedPath.Add(targetNode);

            retracedPath.Reverse();
            return retracedPath;
        }

        private float[,] traversableMap(float[,] oldMap)
        {
            int size = 250;
            float sizex = (float)size / oldMap.GetLength(0);
            float sizez = (float)size / oldMap.GetLength(1);
            float factor = (float)oldMap.GetLength(1) / (float)oldMap.GetLength(0);
            float[,] newMap = new float[size, size];
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    newMap[i, j] = oldMap[(int)(Mathf.Floor(i / (sizex))), (int)(Mathf.Floor(j / sizez))];
                }

            }
            return newMap;
        }

        private Vector3 CheckIfWaypointcollidesAndAdjust(Vector3 wp, float[,] newMap, TerrainInfo tinfo, float[] carsize)
        {
            float[] size = tinfo.getSize(); 
            float xstep = (size[0]) / ((float)newMap.GetLength(0)); //size of step in x direction
            float zstep = (size[1]) / ((float)newMap.GetLength(1)); //size of step in z direction
            double xfactor = ((carsize[0] / xstep)); //number of steps in xdirection corresponding to carsize
            double zfactor = ((carsize[1] / zstep)); //number of steps in zdiretion corresponding to carsize
            Vector3 newVect = wp;
            float[] sqrtcarsize = new float[] {(float)Math.Sqrt((double)carsize[0]), (float)Math.Sqrt((double)carsize[1]) };


            int[] index = new int[] { (int)((wp.x - size[2]) / xstep), (int)((wp.z - size[3]) / zstep) };
            //Debug.Log(index[0] + " " + index[1]);
            //Debug.Log(((int)xfactor + 1));
            /*
            for (int i = 1; i < ((int)xfactor + 1); i++)
            {    
                if (newMap[index[0]+i, index[1]+i] > 0.7)
                {
                    newVect = wp + (new Vector3(-xstep*(float)((xfactor-i+1)), 0, -zstep*(float)((xfactor - i + 1))));
                    break;
                }
                else if (newMap[index[0] - i, index[1] + i] > 0.7)
                {
                    newVect = wp + new Vector3(xstep * (float)((xfactor - i + 1)), 0, -zstep * (float)((xfactor - i + 1)));
                    break;
                }
                else if (newMap[index[0] - i, index[1] - i] > 0.7)
                {
                    newVect = wp + new Vector3(xstep * (float)((xfactor - i + 1)), 0, zstep * (float)((xfactor - i + 1)));
                    break;
                }
                else if (newMap[index[0] + i, index[1] - i] > 0.7)
                {
                    newVect = wp + new Vector3(-xstep * (float)((xfactor - i + 1)), 0, zstep * (float)((xfactor - i + 1)));
                    break;
                }
            }
            */

            return newVect;
        }

        private float[,] InflateObstacles(float[,] newMap, float[] carsize, TerrainInfo tinfo)
        {
            float[] size = tinfo.getSize();
            float xstep = (size[0]) / ((float)newMap.GetLength(0)); //size of step in x direction
            float zstep = (size[1]) / ((float)newMap.GetLength(1)); //size of step in z direction
            double xfactor = ((carsize[0] / xstep)); //number of steps in xdirection corresponding to carsize
            double zfactor = ((carsize[1] / zstep)); //number of steps in zdiretion corresponding to carsize
            float[] sqrtcarsize = new float[] { (float)Math.Sqrt((double)carsize[0]), (float)Math.Sqrt((double)carsize[1]) };


            //Debug.Log(((int)xfactor + 1));
            for (int k = 1; k < newMap.GetLength(0); k++)
            {
                for (int j = 1; j < newMap.GetLength(1); j++)
                {
                    
                    for (int i = 1; i < ((int)xfactor+1); i++)
                    {
                        if (k + i >= newMap.GetLength(0) - 1 || j + i >= newMap.GetLength(1) - 1 || j + i < 0 || j - i < 0 || newMap[k, j] > 0.8f)
                        {
                            break;
                        }
                        else
                        {
                            if (newMap[k + i, j + i] > 0.8)
                            {
                                newMap[k, j] = 0.6f;
                                break;
                            }
                            else if (newMap[k - i, j + i] > 0.8)
                            {
                                newMap[k, j] = 0.6f;
                                break;
                            }
                            else if (newMap[k - i, j - i] > 0.8)
                            {
                                newMap[k, j] = 0.6f;
                                break;
                            }
                            else if (newMap[k + i, j - i] > 0.8)
                            {
                                newMap[k, j] = 0.6f;
                                break;
                            }
                        }
                    }
                }
            }

            return newMap;
        }

        private int large_get_i_index(float x, float[,] newMap, TerrainInfo tinfo)
        {
            float[] size = tinfo.getSize();
            float xstep = (size[0]) / ((float)newMap.GetLength(0)); //size of step in x direction
            float zstep = (size[1]) / ((float)newMap.GetLength(1)); //size of step in z direction
            int index = (int)((x - size[2]) / xstep);
            return index;
        }

        private int large_get_j_index(float z, float[,] newMap, TerrainInfo tinfo)
        {
            float[] size = tinfo.getSize();
            float xstep = (size[0]) / ((float)newMap.GetLength(0)); //size of step in x direction
            float zstep = (size[1]) / ((float)newMap.GetLength(1)); //size of step in z direction
            int index = (int)((z - size[3]) / zstep);
            return index;
        }

        private float large_get_i_pos(int x, float[,] newMap, TerrainInfo tinfo)
        {
            float[] size = tinfo.getSize();
            float xstep = (size[0]) / ((float)newMap.GetLength(0)); //size of step in x direction
            float zstep = (size[1]) / ((float)newMap.GetLength(1)); //size of step in z direction
            float pos = (float)(x * xstep + size[2]);
            return pos;
        }

        private float large_get_j_pos(int z, float[,] newMap, TerrainInfo tinfo)
        {
            float[] size = tinfo.getSize();
            float xstep = (size[0]) / ((float)newMap.GetLength(0)); //size of step in x direction
            float zstep = (size[1]) / ((float)newMap.GetLength(1)); //size of step in z direction
            float pos = (float)(z * xstep + size[3]);
            return pos;
        }

        //small change for improvements

        int wp = 1;
        float rotationangle = 0;
        int towards = 1;
        float lastDistance = 10000000000;
        Vector3 oldPos = Vector3.zero;
        float oldRotAngle = 0;
        float oldRotAngle2 = 0;

        private void FixedUpdate()
        {
            Vector3 carPos = m_Car.transform.position;
            Vector3 velocity = (oldPos - carPos)/ Time.fixedDeltaTime;
            oldPos = carPos;
            Quaternion carRot = m_Car.transform.rotation;
            // Execute your path here
            // ...

            Vector3 axis = new Vector3(1f,0f,0f);
            carRot.ToAngleAxis(out rotationangle, out axis);
            //if(oldRotAngle + (oldRotAngle - oldRotAngle2) < (float)(Math.PI))
            if(carRot.y <= 0)
            {
                rotationangle = rotationangle / 360 * 2 * (float)(Math.PI);
                oldRotAngle = rotationangle;
                oldRotAngle2 = oldRotAngle;

            }
            else
            {
                rotationangle = 2 * (float)(Math.PI) - rotationangle / 360 * 2 * (float)(Math.PI);
                oldRotAngle = rotationangle;
                oldRotAngle2 = oldRotAngle;
            }
            
            //Debug.Log(rotationangle + "rotationangle");
            //Debug.Log(carRot.w + "quatval  ");
            //Debug.Log(carRot.x + "quatval  ");
            //Debug.Log(carRot.y + "quatval  ");
            //Debug.Log(carRot.z + "quatval  ");




            // this is how you access information about the terrain from the map
            int i = terrain_manager.myInfo.get_i_index(transform.position.x);
            int j = terrain_manager.myInfo.get_j_index(transform.position.z);
            float grid_center_x = terrain_manager.myInfo.get_x_pos(i);
            float grid_center_z = terrain_manager.myInfo.get_z_pos(j);

            Debug.DrawLine(transform.position, new Vector3(grid_center_x, 0f, grid_center_z));

            // this is how you access information about the terrain from a simulated laser range finder
            RaycastHit hit;
            float maxRange = 50f;
            if (Physics.Raycast(transform.position + transform.up, transform.TransformDirection(Vector3.forward), out hit, maxRange))
            {
                Vector3 closestObstacleInFront = transform.TransformDirection(Vector3.forward) * hit.distance;
                Debug.DrawRay(transform.position, closestObstacleInFront, Color.yellow);
                //Debug.Log("Did Hit");
            }

            float distWP = PathSegment.getDistToWP(my_path[wp], carPos.z, carPos.x);
            float distRay = PathSegment.getDistToRay(pathAngles[wp], my_path[wp], carPos.z, carPos.x);
            float biasRay = 0;
            if (distRay > 0.3f)
            {
                //biasRay = 20;
            }

            if (pathAngles[wp] == 100f)
            {
                float angleGoal = PathSegment.getAngleCarGoal(carPos, my_path[wp]);
                distRay = PathSegment.getDistToRay(angleGoal, my_path[wp], carPos.z, carPos.x);
                float angleDiff = -rotationangle - angleGoal;
                //Debug.Log(rotationangle + " rotationangle");
                //Debug.Log(angleGoal);
                //Debug.Log(distWP);
                float turningratio = 10*angleDiff / distWP;
                m_Car.Move(turningratio, 0.8f, 0f, 0f);
            }
            else
            {

                if (distRay > lastDistance)
                {
                    towards = -1;
                    distRay = 0.4f;
                    biasRay = 20;
                }
                lastDistance = distRay;



                //Debug.Log(distRay + " 0");
                //Debug.Log(lastDistance);


                float dir = 1;
                //float angleWP = -PathSegment.getAngleForMiddle(carPos, my_path[wp]);
                if (Math.Abs((double)(rotationangle - pathAngles[wp])) > Math.PI)
                {
                    dir = -1;
                }
                float angleDiff = (1 + biasRay)*(rotationangle - pathAngles[wp]);
                Debug.Log(pathAngles[wp]);
                //Debug.Log(angleDiff + " anglediff");
                Debug.Log(rotationangle + " rotationangle");
                //Debug.Log(distRay);


                float turningratio = (angleDiff*dir) / (float)(Math.Pow((double)distRay, 2));
                if (PathSegment.getDistToWP(my_path[wp], carPos.z, carPos.x) < 1)
                {
                    wp++;
                    towards = 1;
                    lastDistance = 100000;
                }



                //Debug.Log(pathAngles[wp] + " path");
                //Debug.Log(pathAngles[0] + " path");
                //Debug.Log(pathAngles[2] + " path");
                //Debug.Log(pathAngles[3] + " path");

                //Debug.Log(rotationangle);
                //Debug.Log(wp);


                // this is how you control the car
                float acc = 0.8f;
                if (velocity.magnitude > 5f)
                {
                    acc = 0f;
                }
                m_Car.Move(turningratio, acc, 0f, 0f);
            }

        }
    }


    public class PathSegment
    {
        public static float[] circle1 = new float[] { 1, 1, 1, 0, 90 }; //x, z, radius, start angle, end angle
        public static float[] circle2 = new float[] { 1, 1, 1, 0, 90 };
        public static float[] line = new float[] { 1, 1, 2, 2 }; //x1,z1,x2,z2

        public static void createPathSeg(float startAngle, float endAngle, Vector3 wp1, Vector3 wp2, float cutoffAngle)
        {
            int r = 10;
            circle1 = new float[] { wp1.x - (float)(r*Math.Cos((double)startAngle)), wp1.z - (float)(r * Math.Sin((double)startAngle)), r, startAngle, cutoffAngle };
            circle2 = new float[] { wp2.x - (float)(r * Math.Cos((double)(endAngle))), wp2.z - (float)(r * Math.Sin((double)(endAngle))), r, cutoffAngle, endAngle };
            //Debug.Log(circle2[0] + " " + circle2[1] + " " + wp2.x + " " + wp2.z);
            line = new float[] { wp1.x - (float)(r * Math.Cos((double)startAngle)) + (float)(r * Math.Cos((double)circle1[4])),
                wp1.z - (float)(r * Math.Sin((double)startAngle)) + (float)(r * Math.Sin((double)circle1[4])),
            wp2.x - (float)(r * Math.Cos((double)(circle2[4]))) + (float)(r * Math.Cos((double)circle2[3])),
            wp2.z - (float)(r * Math.Sin((double)(circle2[4]))) + (float)(r * Math.Sin((double)circle2[3]))};

            //line
            Debug.DrawLine(new Vector3(line[0], 0, line[1]), new Vector3(line[2], 0, line[3]), Color.red, 100f);

            //circle1
            Debug.DrawLine(new Vector3(circle1[0] + circle1[2]*(float)Math.Cos((double)circle1[3]), 0, circle1[1] + circle1[2] * (float)Math.Sin((double)circle1[3])),
                new Vector3(circle1[0] + circle1[2] * (float)Math.Cos((double)circle1[4]), 0, circle1[1] + circle1[2] * (float)Math.Sin((double)circle1[4])), Color.blue, 100f);

            Debug.DrawLine(new Vector3(circle1[0] + circle1[2] * (float)Math.Cos((double)circle1[3]), 0, circle1[1] + circle1[2] * (float)Math.Sin((double)circle1[3])),
                new Vector3(circle1[0], 0, circle1[1]), Color.blue, 100f);

            Debug.DrawLine(new Vector3(circle1[0], 0, circle1[1]),
                new Vector3(circle1[0] + circle1[2] * (float)Math.Cos((double)circle1[4]), 0, circle1[1] + circle1[2] * (float)Math.Sin((double)circle1[4])), Color.blue, 100f);

            //circle2
            Debug.DrawLine(new Vector3(circle2[0] + circle2[2] * (float)Math.Cos((double)circle2[3]), 0, circle2[1] + circle2[2] * (float)Math.Sin((double)circle2[3])),
                new Vector3(circle2[0] + circle2[2] * (float)Math.Cos((double)circle2[4]), 0, circle2[1] + circle2[2] * (float)Math.Sin((double)circle2[4])), Color.blue, 100f);

            Debug.DrawLine(new Vector3(circle2[0], 0, circle2[1]),
                new Vector3(circle2[0] + circle2[2] * (float)Math.Cos((double)circle2[4]), 0, circle2[1] + circle2[2] * (float)Math.Sin((double)circle2[4])), Color.blue, 100f);

            Debug.DrawLine(new Vector3(circle2[0] + circle2[2] * (float)Math.Cos((double)circle2[3]), 0, circle2[1] + circle2[2] * (float)Math.Sin((double)circle2[3])),
                new Vector3(circle2[0], 0, circle2[1]), Color.blue, 100f);

        }

        public static void createVertexLine(float angle, Vector3 wp1)
        {
            Debug.DrawLine(wp1 + new Vector3((float)(100 * Math.Cos(angle)), 0, (float)(100 * Math.Sin(angle))), wp1 + new Vector3(-(float)(100 * Math.Cos(angle)), 0, -(float)(100 * Math.Sin(angle))), Color.red, 100f);
        }

        public static float getDistToRay(float angle, Vector3 wp1, float carz, float carx)
        {
            angle = angle - (float)(Math.PI / 2);
            return (float)(Math.Abs((double)(Math.Cos((double)angle) * (wp1.z - carz) - Math.Sin((double)angle) * (wp1.x - carx))));
        }

        public static float getDistToWP(Vector3 wp1, float carz, float carx)
        {
            return (float)(Math.Sqrt(Math.Pow((double)(carx-wp1.x),2) + Math.Pow((double)(carz - wp1.z), 2)));
        }

        public static float getAngleForMiddle(Vector3 wp1, Vector3 wp2)
        {
            float val = (float)(Math.Atan(((double)((wp1.x - wp2.x) / (wp2.z - wp1.z)))));// 90 degrees defined as u
            //Debug.Log(val + ":angle");
            //Debug.Log(wp1.x - wp2.x);
            //Debug.Log(wp1.z - wp2.z);
            if(wp1.x - wp2.x == 0)
            {
                if (wp2.z - wp1.z < 0)
                {
                    return (float)Math.PI;
                }else
                {
                    return 0f;
                }
            }else if((wp1.x - wp2.x)  < 0)
            {
                if ((wp2.z - wp1.z) < 0)
                {
                    val += (float)Math.PI;
                    return val;
                }
                val += 2*(float)Math.PI;
                return val;
            }

            if (val < 0)
            {
                val += (float)Math.PI;
            }
            return val;
        }

        public static float getAngleCarGoal(Vector3 wp1, Vector3 wp2)
        {
            float val = (float)(Math.Atan(((double)((wp1.z - wp2.z) / (wp1.x - wp2.x)))));// 90 degrees defined as u
            Debug.Log(val + ":angle");
            Debug.Log(wp1.x - wp2.x);
            Debug.Log(wp1.z - wp2.z);
            if (val < 0)
            {
                //val += 2*(float)Math.PI;
            }
            return val;
        }
    }

    public class PathPoint
    {
        public PathPoint parentNode;
        //Indices
        public int i;
        public int j;

        //Distance to Goal
        public float h;
        //Distance travelled
        public float g;
        //Total Fitness
        public float f;

        //direction
        public Vector2 direction;

        public void Start()
        {
            i = 0;
            j = 0;
            h = 0.0f;
            g = 0.0f;
            f = 0.0f;
            direction = new Vector2(0f, 0f);
            parentNode = this;
        }

        public void SetParent(PathPoint parent)
        {
            this.parentNode = parent;
        }
        public PathPoint GetParent()
        {
            return this.parentNode;
        }
    }
}

