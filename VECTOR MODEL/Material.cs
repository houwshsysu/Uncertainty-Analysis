using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ESRI.ArcGIS.Geometry;
using System.Reflection;
using ESRI.ArcGIS.Geodatabase;

namespace VECTOR_MODEL
{
    public class Material
    {
        private static object _missing = Type.Missing;

        public static void Clone(IVector3D ObjectVector, IVector3D GetVector)
        {
            GetVector.XComponent = ObjectVector.XComponent;
            GetVector.YComponent = ObjectVector.YComponent;
            GetVector.ZComponent = ObjectVector.ZComponent;
        }

        //根据两点构建3D向量
        public static IVector3D CreateVector3DTwoPoints(IPoint toPoint, IPoint fromPoint)
        {
            if (fromPoint.Z.Equals(Double.NaN))
            {
                fromPoint.Z = 0;
            }
            if (toPoint.Z.Equals(Double.NaN))
            {
                toPoint.Z = 0;
            }
            IVector3D vector3D = new Vector3DClass();
            vector3D.ConstructDifference(toPoint, fromPoint);
            return vector3D;
        }

        public static class GeometryUtilities
        {
            private static object _missing = Type.Missing;

            public static void MakeZAware(IGeometry geometry)
            {
                IZAware zAware = geometry as IZAware;
                zAware.ZAware = true;
            }

            public static IVector3D ConstructVector3D(double xComponent, double yComponent, double zComponent)
            {
                IVector3D vector3D = new Vector3DClass();
                vector3D.SetComponents(xComponent, yComponent, zComponent);

                return vector3D;
            }

            public static double GetRadians(double decimalDegrees)
            {
                return decimalDegrees * (Math.PI / 180);
            }

            public static IPoint ConstructPoint3D(double x, double y, double z)
            {
                IPoint point = ConstructPoint2D(x, y);
                point.Z = z;
                MakeZAware(point as IGeometry);
                return point;
            }

            public static IPoint ConstructPoint2D(double x, double y)
            {
                IPoint point = new PointClass();
                point.X = x;
                point.Y = y;

                return point;
            }
        }

        //曲率熵
        //参数已排序                 
        public static double Entropy(ITinNode ENode, double[] ECurvature,ITinAdvanced2 Etin)
        {
            //曲率值绝对值化
            int ELength = ECurvature.Length;
            double[] AbsE = new double[ELength];
            for (int i = 0; i < ELength; i++)
            {
                AbsE[i] = Math.Abs(ECurvature[i]);
            }
            //得到该点某曲率绝对值
            int u = MNodeLocation(ENode);
            double Ki = AbsE[u];
            //得到一阶领域点曲率绝对值
            ITinNodeArray ENodeArray = GetIncidentNodes(ENode,VectorModel.NodeRange);
            int ENodeCount = ENodeArray.Count;
            double[] Kj = new double[ENodeCount];
            for (int i = 0; i < ENodeCount; i++)
            {
                ITinNode ENearNode = ENodeArray.get_Element(i);
                int v = MNodeLocation(ENearNode);
                Kj[i] = AbsE[v];
            }
            double KSum = Ki;
            for (int i = 0; i < ENodeCount; i++)
            {
                KSum = KSum + Kj[i];
            }
            double Pi = Ki / KSum;
            double[] Pj = new double[ENodeCount];
            for (int i = 0; i < ENodeCount; i++)
            {
                Pj[i] = Kj[i] / KSum;
            }
            double Hi = -Pi * Math.Log(Pi, 2);
            for (int i = 0; i < ENodeCount; i++)
            {
                Hi = Hi - Pj[i] * Math.Log(Pj[i], 2);
            }
            return Hi;
        }

        //得到箭头
        public static IGeometry GetArrow()
        {
            const double ConeBaseDegrees = 360.0;
            const int ConeBaseDivisions = 36;
            const double VectorComponentOffset = 0.0000001;
            const double ConeBaseZ = 0.0;
            double ConeApexZ = 0;
            double ConeBaseRadius = 0;

            ConeBaseRadius = 20;                             //范围                           //正常
            ConeApexZ = 100;                                   //尖度  

            //Vector3D: Cone, TriangleFan With 36 Vertices

            IGeometryCollection multiPatchGeometryCollection = new MultiPatchClass();

            IPointCollection triangleFanPointCollection = new TriangleFanClass();

            //Set Cone Apex To (0, 0, ConeApexZ)

            IPoint coneApexPoint = GeometryUtilities.ConstructPoint3D(0, 0, ConeApexZ);

            //Add Cone Apex To Triangle Fan

            triangleFanPointCollection.AddPoint(coneApexPoint, ref _missing, ref _missing);

            //Define Upper Portion Of Axis Around Which Vector Should Be Rotated To Generate Cone Base Vertices

            IVector3D upperAxisVector3D = GeometryUtilities.ConstructVector3D(0, 0, 10);

            //Define Lower Portion of Axis Around Which Vector Should Be Rotated To Generate Cone Base Vertices

            IVector3D lowerAxisVector3D = GeometryUtilities.ConstructVector3D(0, 0, -10);

            //Add A Slight Offset To X or Y Component Of One Of Axis Vectors So Cross Product Does Not Return A Zero-Length Vector

            lowerAxisVector3D.XComponent += VectorComponentOffset;

            //Obtain Cross Product Of Upper And Lower Axis Vectors To Obtain Normal Vector To Axis Of Rotation To Generate Cone Base Vertices

            IVector3D normalVector3D = upperAxisVector3D.CrossProduct(lowerAxisVector3D) as IVector3D;

            //Set Normal Vector Magnitude Equal To Radius Of Cone Base

            normalVector3D.Magnitude = ConeBaseRadius;

            //Obtain Angle Of Rotation In Radians As Function Of Number Of Divisions Within 360 Degree Sweep Of Cone Base

            double rotationAngleInRadians = GeometryUtilities.GetRadians(ConeBaseDegrees / ConeBaseDivisions);

            for (int i = 0; i < ConeBaseDivisions; i++)
            {
                //Rotate Normal Vector Specified Rotation Angle In Radians Around Either Upper Or Lower Axis

                normalVector3D.Rotate(-1 * rotationAngleInRadians, upperAxisVector3D);

                //Construct Cone Base Vertex Whose XY Coordinates Are The Sum Of Apex XY Coordinates And Normal Vector XY Components

                IPoint vertexPoint = GeometryUtilities.ConstructPoint3D(coneApexPoint.X + normalVector3D.XComponent,
                                                                      coneApexPoint.Y + normalVector3D.YComponent,
                                                                      ConeBaseZ);

                //Add Vertex To TriangleFan

                triangleFanPointCollection.AddPoint(vertexPoint, ref _missing, ref _missing);
            }

            //Re-Add The Second Point Of The Triangle Fan (First Vertex Added) To Close The Fan

            triangleFanPointCollection.AddPoint(triangleFanPointCollection.get_Point(1), ref _missing, ref _missing);

            //Add TriangleFan To MultiPatch

            multiPatchGeometryCollection.AddGeometry(triangleFanPointCollection as IGeometry, ref _missing, ref _missing);

            return multiPatchGeometryCollection as IGeometry;
        }

        public static void MakeZAware(IGeometry geometry)
        {
            IZAware zAware = geometry as IZAware;
            zAware.ZAware = true;
        }

        public static IVector3D ConstructVector3D(double xComponent, double yComponent, double zComponent)
        {
            IVector3D vector3D = new Vector3DClass();
            vector3D.SetComponents(xComponent, yComponent, zComponent);

            return vector3D;
        }

        public static double GetRadians(double decimalDegrees)
        {
            return decimalDegrees * (Math.PI / 180);
        }

        public static IPoint ConstructPoint3D(double x, double y, double z)
        {
            IPoint point = ConstructPoint2D(x, y);
            point.Z = z;

            MakeZAware(point as IGeometry);

            return point;
        }

        public static IPoint ConstructPoint2D(double x, double y)
        {
            IPoint point = new PointClass();
            point.X = x;
            point.Y = y;

            return point;
        }

        #region 排序函数
        public static int MNodeLocation(ITinNode NodeL)
        {
            int v = (int)((NodeL.X - (GetOriginal.XMin - (GetOriginal.RasterSize / 2))) / GetOriginal.RasterSize) + GetOriginal.RowCount * (int)((NodeL.Y - (GetOriginal.YMin - (GetOriginal.RasterSize / 2))) / (GetOriginal.RasterSize));
            return v;
        }

        public static int MPointLocation(IPoint PointL)
        {
            int v = (int)((PointL.X - (GetOriginal.XMin - (GetOriginal.RasterSize / 2))) / GetOriginal.RasterSize) + GetOriginal.RowCount * (int)((PointL.Y - (GetOriginal.YMin - (GetOriginal.RasterSize / 2))) / (GetOriginal.RasterSize));
            return v;
        }
        #endregion

        #region 备用邻域函数
        //////最多只能二阶
        //备用邻域函数
        //节点
        public static ITinNodeArray GetIncidentNodes(ITinNode Node, double range)
        {
            ITinNodeArray Incident1 = Node.GetAdjacentNodes();
            int IncidentNodeCount1 = Incident1.Count;
            ITinNodeArray[] Incident2 = new ITinNodeArray[IncidentNodeCount1];
            //设置的长度100
            ITinNode[] IsRepeat = new ITinNode[1000];
            int indec = 0;
            for (int i = 0; i < IncidentNodeCount1; i++)
            {
                ITinNode node1 = Incident1.get_Element(i);
                ITinNodeArray Node2 = node1.GetAdjacentNodes();
                Incident2[i] = Node2;
            }
            for (int i = 0; i < IncidentNodeCount1; i++)
            {
                int IncidentNodeCount2 = Incident2[i].Count;
                for (int j = 0; j < IncidentNodeCount2; j++)
                {
                    IsRepeat[indec] = Incident2[i].get_Element(j);
                    indec++;
                }
            }
            for (int i = 0; i < 1000; i++)
            {
                if (IsRepeat[i] == null)
                {
                    break;
                }
                bool IsNeed = true;
                for (int j = 0; j < Incident1.Count; j++)
                {
                    if ((Incident1.get_Element(j).Index == IsRepeat[i].Index) || (Node.Index == IsRepeat[i].Index) || Math.Abs(Node.X - IsRepeat[i].X) > (GetOriginal.RasterSize * range) || (Math.Abs(Node.Y - IsRepeat[i].Y) > (GetOriginal.RasterSize * range)))
                    {
                        IsNeed = false;
                        break;
                    }
                }
                if (IsNeed == true)
                {
                    Incident1.Add(IsRepeat[i]);
                }
            }

            //排除四个角的4点
            int FinCount = Incident1.Count;
            for (int i = 0; i < FinCount; i++)
            {
                ITinNode CheckNode = Incident1.get_Element(i);
                int CheckID = CheckNode.Index;
                if (CheckID < 5)
                {
                    Incident1.Remove(i);
                    FinCount--;
                    i--;
                }
            }
            return Incident1;
        }

        //三角形
        public static ITinTriangleArray GetIncidentTriangles(ITinNode Node, double range)
        {
            ITinNodeArray IncidentNode1 = Node.GetAdjacentNodes();
            ITinTriangleArray IncidentTriangle1 = Node.GetIncidentTriangles();
            int IncidentNodeCount1 = IncidentNode1.Count;
            ITinTriangleArray[] IncidentTriangle2 = new ITinTriangleArray[IncidentNodeCount1];
            //设置的长度1000
            ITinTriangle[] IsRepeat = new ITinTriangle[1000];
            int indec = 0;
            for (int i = 0; i < IncidentNodeCount1; i++)
            {
                ITinNode node1 = IncidentNode1.get_Element(i);
                ITinTriangleArray Triangle2 = node1.GetIncidentTriangles();
                IncidentTriangle2[i] = Triangle2;
            }
            for (int i = 0; i < IncidentNodeCount1; i++)
            {
                int IncidentNodeCount2 = IncidentTriangle2[i].Count;
                for (int j = 0; j < IncidentNodeCount2; j++)
                {
                    IsRepeat[indec] = IncidentTriangle2[i].get_Element(j);
                    indec++;
                }
            }
            for (int i = 0; i < 1000; i++)
            {
                if (IsRepeat[i] == null)
                {
                    break;
                }
                IPoint Cpoint = new PointClass();
                IsRepeat[i].QueryCentroid(Cpoint);
                bool IsNeed = true;
                for (int j = 0; j < IncidentTriangle1.Count; j++)
                {
                    if ((IsRepeat[i].Index == IncidentTriangle1.get_Element(j).Index) || Math.Abs(Cpoint.X - Node.X) > (GetOriginal.RasterSize * range) || Math.Abs(Cpoint.Y - Node.Y) > (GetOriginal.RasterSize * range))
                    {
                        IsNeed = false;
                        break;
                    }
                }
                if (IsNeed == true)
                {
                    IncidentTriangle1.Add(IsRepeat[i]);
                }
            }
            return IncidentTriangle1;
        }
        #endregion
    }
}

