
/* 得到单个曲面上的矢量参量表征体系
 */

using System;
using System.Collections;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using ESRI.ArcGIS.Analyst3D;
using ESRI.ArcGIS.Geodatabase;
using ESRI.ArcGIS.Carto;
using ESRI.ArcGIS.DataSourcesFile;
using ESRI.ArcGIS.esriSystem;
using ESRI.ArcGIS.Geometry;
using System.IO;

namespace VECTOR_MODEL
{
    class GetOriginal
    {
        //---------------------------------------------------------------------------------------------
        //工作空间路径
        private string workSpace = @"C:\Users\B410\Desktop\beforeSchool\optimization\FdOrigin3-1-1.5.gdb";
        //曲面路径
        private string tinPath = @"C:\Users\B410\Desktop\C&G\Data\fd";
        //txt文件路径
        private string pointsFile = @"C:\Users\B410\Desktop\beforeSchool\BBBB\FdSinglePointsFile3-1-1.5.txt";
        //txt文件路径
        private string exposepath = @"C:\Users\B410\Desktop\C&G\Data\briefBore\expose.txt";
        //栅格单元边长 
        public static double RasterSize = 5;
        //---------------------------------------------------------------------------------------------
        public static long TriangleCount;
        public static long NodeCount;
        IVector3D[] NormalVectorOfNode;
        IVector3D[] NormalVectorOfSurface;
        IPoint[] Centroid;
        ITinAdvanced2 tin;
        IWorkspace Workspace;
        ITinNode[] Node;
        IPoint[] Point;
        double[] KG;
        double[] H;
        double[] k1;
        double[] k2;
        IVector3D[] e1;
        IVector3D[] e2;
        double[] HEntropy;
        double[] GEntropy;
        //曲率正负判定
        int IsRegular;
        //完成排序后的参数
        public static IVector3D[] ONV;
        public static IVector3D[] Oe1;
        public static IVector3D[] Oe2;
        public static double[] Ok1;
        public static double[] Ok2;
        public static double[] OG;
        public static double[] OM;
        public static double[] OME;
        public static double[] OGE;
        public static double[] OSTD;
        public static IPoint[] OPoint;
        public static ITinNode[] ONode;
        //曲面信息
        public static int RowCount;
        public static double XMin;
        public static double XMax;
        public static double YMin;
        public static double YMax;
        //要素类声明
        IFields NodeFields;
        IFeatureClass NodeFeatureClass;
        IFields UnitNormalVectorOfNodeLinefields;
        IFields UnitNormalVectorOfNodeArrowfields;
        IFeatureClass UnitNormalVectorOfNodeLineFeatureClass;
        IFeatureClass UnitNormalVectorOfNodeArrowFeatureClass;
        IFields NodeVectorLinefieldsOne;
        IFields NodeVectorArrowfieldsOne;
        IFeatureClass NodeVectorLineFeatureClassOne;
        IFeatureClass NodeVectorArrowFeatureClassOne;
        IFields NodeVectorLinefieldsTwo;
        IFields NodeVectorArrowfieldsTwo;
        IFeatureClass NodeVectorLineFeatureClassTwo;
        IFeatureClass NodeVectorArrowFeatureClassTwo;

        int[] ChooseBore;

        private static object _missing = Type.Missing;

        // 获得TIN数据集
        public void GetTinData()
        {
            Type factoryType = Type.GetTypeFromProgID("esriDataSourcesGDB.FileGDBWorkspaceFactory");
            IWorkspaceFactory workspaceFactory = (IWorkspaceFactory)Activator.CreateInstance(factoryType);
            Workspace = workspaceFactory.OpenFromFile(workSpace, 0);
            // Open an existing TIN, edit the following path as appropriate.
            tin = new TinClass();
            tin.Init(tinPath);
            //得到和tin网直接相关的数据
            TriangleCount = tin.DataTriangleCount;
            NodeCount = tin.NodeCount;
            Node = new ITinNode[NodeCount];
            Point = new IPoint[NodeCount];
            IEnvelope Envelope = tin.Extent;
            XMin = Envelope.XMin;
            XMax = Envelope.XMax;
            YMin = Envelope.YMin;
            YMax = Envelope.YMax;
            //O数组
            OPoint = new IPoint[NodeCount];
            ONode = new ITinNode[NodeCount];
            ONV = new IVector3D[NodeCount];
            Oe1 = new IVector3D[NodeCount];
            Oe2 = new IVector3D[NodeCount];
            Ok1 = new double[NodeCount];
            Ok2 = new double[NodeCount];
            OG = new double[NodeCount];
            OM = new double[NodeCount];
            OME = new double[NodeCount];
            OGE = new double[NodeCount];
            OSTD = new double[NodeCount];

            //列数
            RowCount = (int)((XMax - XMin + 0.5 * RasterSize) / RasterSize) + 1;
            for (int i = 4; i < NodeCount; i++)
            {
                //所有点集
                Node[i] = tin.GetNode(i + 1);
                Point[i] = new PointClass();
                Node[i].QueryAsPoint(Point[i]);
                Material.MakeZAware(Point[i] as IGeometry);
                int u = PointLocation(Point[i]);
                OPoint[u] = Point[i];
                Material.MakeZAware(OPoint[u]);
                ONode[u] = Node[i];
            }
        }

        #region NormalVector
        //获取每个三角形的法矢量
        public void GetSurfaceNormalVector()
        {
            ITinAdvanced TinAdvance = tin as ITinAdvanced;
            Centroid = new IPoint[TriangleCount];
            NormalVectorOfSurface = new IVector3D[TriangleCount];
            for (int i = 0; i < TriangleCount; i++)
            {
                Centroid[i] = new PointClass();
                Material.MakeZAware(Centroid[i] as IGeometry);
                NormalVectorOfSurface[i] = new Vector3DClass();
                ITinTriangle TinTriangle = TinAdvance.GetTriangle(i + 1);
                if (TinTriangle.IsEmpty)
                    continue;
                //得到三角形形心
                TinTriangle.QueryCentroid(Centroid[i]);
                //得到基于面的单位法矢量
                TinTriangle.QueryNormal(NormalVectorOfSurface[i]);
                if (NormalVectorOfSurface[i].IsEmpty)
                    continue;
                NormalVectorOfSurface[i].Normalize();
            }
        }

        //获取网格节点p的法矢量
        public void GetNodeNormalVector()
        {
            ITinAdvanced TinAdvance = tin as ITinAdvanced;
            NormalVectorOfNode = new IVector3D[NodeCount];
            //对每个节点进行循环
            for (int i = 4; i < NodeCount; i++)
            {
                NormalVectorOfNode[i] = new Vector3DClass();
                ITinTriangleArray IncidentTriangles = Material.GetIncidentTriangles(Node[i], VectorModel.CentriodRange);
                int IncidentTriangleCount = IncidentTriangles.Count;
                //对p节点的一阶邻域节点循环
                for (int j = 0; j < IncidentTriangleCount; j++)
                {
                    IPoint IncidentCentroid = new PointClass();
                    Material.MakeZAware(IncidentCentroid as IGeometry);
                    IVector3D IncidentUnitNormalVectorOfSurface = new Vector3DClass();
                    ITinTriangle IncidentTinTriangle = IncidentTriangles.get_Element(j);
                    if (IncidentTinTriangle.IsEmpty)
                        continue;
                    IncidentTinTriangle.QueryNormal(IncidentUnitNormalVectorOfSurface);
                    if (IncidentUnitNormalVectorOfSurface.IsEmpty)
                        continue;
                    //得到p节点相邻三角形的形心和单位法矢量,得到节点到形心的距离
                    IncidentTinTriangle.QueryCentroid(IncidentCentroid);
                    double NCdistance = System.Math.Sqrt((IncidentCentroid.X - Node[i].X) * (IncidentCentroid.X - Node[i].X) +
                    (IncidentCentroid.Y - Node[i].Y) * (IncidentCentroid.Y - Node[i].Y) + (IncidentCentroid.Z - Node[i].Z) * (IncidentCentroid.Z - Node[i].Z));
                    //得到节点的法矢量分量,合成法矢量
                    IncidentUnitNormalVectorOfSurface.Magnitude = (1 / NCdistance);
                    if ((NormalVectorOfNode[i] == null) || (NormalVectorOfNode[i].IsEmpty))
                    {
                        NormalVectorOfNode[i] = IncidentUnitNormalVectorOfSurface;
                    }
                    else
                    {
                        NormalVectorOfNode[i].ConstructAddVector(NormalVectorOfNode[i], IncidentUnitNormalVectorOfSurface);
                    }
                }
                if (NormalVectorOfNode[i].IsEmpty)
                    continue;
                NormalVectorOfNode[i].Normalize();
                //O数组
                int u = (int)((Point[i].X - (XMin - (RasterSize / 2))) / RasterSize) + RowCount * (int)((Point[i].Y - (YMin - (RasterSize / 2))) / (RasterSize));
                ONV[u] = NormalVectorOfNode[i];
            }
        }
        #endregion  NormalVector

        #region DifferentialGeometry
        public void GetNodeCurvature()
        {
            //将向量 pi−p 投影到 p 点所在的切平面,得到切平面的切向单位向量ti
            ITinAdvanced TinAdvance = tin as ITinAdvanced;
            IVector3D[] tid = new IVector3D[NodeCount];
            double[] a = new double[NodeCount];
            double[] b = new double[NodeCount];
            double[] c = new double[NodeCount];
            IVector3D[] E1 = new IVector3D[NodeCount];
            IVector3D[] E2 = new IVector3D[NodeCount];
            double[] a11 = new double[NodeCount];
            double[] a12 = new double[NodeCount];
            double[] a21 = new double[NodeCount];
            double[] a22 = new double[NodeCount];
            double[] a13 = new double[NodeCount];
            double[] a23 = new double[NodeCount];
            double[] O = new double[NodeCount];
            KG = new double[NodeCount];
            H = new double[NodeCount];
            k1 = new double[NodeCount];
            k2 = new double[NodeCount];
            e1 = new IVector3D[NodeCount];
            e2 = new IVector3D[NodeCount];
            HEntropy = new double[NodeCount];
            GEntropy = new double[NodeCount];

            //对每个节点进行循环
            for (int i = 4; i < NodeCount; i++)
            {
                ITinNodeArray AdjacentNodesArray = Material.GetIncidentNodes(Node[i], VectorModel.NodeRange);
                int AdjacentNodesCount = AdjacentNodesArray.Count;
                //沿着edge的切向向量ti
                IVector3D[] ti = new IVector3D[AdjacentNodesCount];
                //沿ti方向的法曲率Knti
                double[] Knti = new double[AdjacentNodesCount];

                //对p节点的一阶邻域循环得到edge曲率
                for (int j = 0; j < AdjacentNodesCount; j++)
                {
                    //求ti
                    ITinNode AdjacentNode = AdjacentNodesArray.get_Element(j);
                    //转换Node为Point
                    IPoint AdjacentPoint = new PointClass();
                    AdjacentNode.QueryAsPoint(AdjacentPoint);
                    Material.MakeZAware(AdjacentPoint as IGeometry);
                    //(pi-p)
                    IVector3D EdgeVectorsOfNode = Material.CreateVector3DTwoPoints(AdjacentPoint, Point[i]);
                    //<pi-p,N>N
                    IVector3D ProjectNormalVector = new Vector3DClass();
                    Material.Clone(NormalVectorOfNode[i], ProjectNormalVector);
                    IVector3D NormalVectorOfNode2 = new Vector3DClass();
                    Material.Clone(NormalVectorOfNode[i], NormalVectorOfNode2);
                    if (ProjectNormalVector.IsEmpty)
                        continue;
                    if (NormalVectorOfNode[i].IsEmpty)
                        continue;
                    ProjectNormalVector.Magnitude = EdgeVectorsOfNode.DotProduct(NormalVectorOfNode2);
                    //ti=(pi-p)-<pi-p,N>N
                    ti[j] = new Vector3DClass();
                    ti[j].ConstructSubtractVector(EdgeVectorsOfNode, ProjectNormalVector);
                    ti[j].Normalize();
                    //Ni-N
                    IVector3D NiN = new Vector3DClass();
                    if (NormalVectorOfNode[(AdjacentNode.Index) - 1] == null)
                        continue;
                    if (NormalVectorOfNode[(AdjacentNode.Index) - 1].IsEmpty)
                        continue;
                    if (NormalVectorOfNode[i].IsEmpty)
                        continue;
                    NiN.ConstructSubtractVector(NormalVectorOfNode[(AdjacentNode.Index) - 1], NormalVectorOfNode[i]);
                    //Kn(ti)
                    Knti[j] = -(EdgeVectorsOfNode.DotProduct(NiN)) / (EdgeVectorsOfNode.DotProduct(EdgeVectorsOfNode));
                }
                //得到最大的Knti曲率Kntid
                double Kntid = Knti[0];
                for (int j = 1; j < AdjacentNodesCount; j++)
                {
                    if (Knti[j] > Kntid)
                    {
                        Kntid = Knti[j];
                    }
                    a[i] = Kntid;
                }
                //得到节点p的最大Kntid所对应的ti，即tid
                for (int j = 0; j < AdjacentNodesCount; j++)
                {
                    if (Knti[j] == Kntid)
                    {
                        tid[i] = ti[j];
                        break;
                    }
                }
                //得到坐标轴e1,e2,夹角cos(&i),sin(&i),以及系数aij
                E1[i] = tid[i];
                if (NormalVectorOfNode[i].IsEmpty)
                    continue;
                E2[i] = new Vector3DClass();
                E2[i].ConstructCrossProduct(E1[i], NormalVectorOfNode[i]);
                E2[i].Normalize();

                //计算各项系数
                double[] A11 = new double[AdjacentNodesCount];
                double[] A12 = new double[AdjacentNodesCount];
                double[] A21 = new double[AdjacentNodesCount];
                double[] A22 = new double[AdjacentNodesCount];
                double[] A13 = new double[AdjacentNodesCount];
                double[] A23 = new double[AdjacentNodesCount];

                for (int j = 0; j < AdjacentNodesCount; j++)
                {
                    double cos = ti[j].DotProduct(E1[i]);
                    double sin = ti[j].DotProduct(E2[i]);
                    A11[j] = cos * cos * sin * sin;
                    A12[j] = cos * sin * sin * sin;
                    A21[j] = A12[j];
                    A22[j] = sin * sin * sin * sin;
                    A13[j] = (Knti[j] - a[i] * cos * cos) * cos * sin;
                    A23[j] = (Knti[j] - a[i] * cos * cos) * sin * sin;
                }

                a11[i] = a12[i] = a21[i] = a22[i] = a13[i] = a23[i] = 0;
                for (int j = 0; j < AdjacentNodesCount; j++)
                {
                    a11[i] += A11[j];
                    a12[i] += A12[j];
                    a21[i] += A21[j];
                    a22[i] += A22[j];
                    a13[i] += A13[j];
                    a23[i] += A23[j];
                }

                //得到系数b，c
                b[i] = (a13[i] * a22[i] - a23[i] * a12[i]) / (a11[i] * a22[i] - a12[i] * a12[i]);
                c[i] = (a11[i] * a23[i] - a12[i] * a13[i]) / (a11[i] * a22[i] - a12[i] * a12[i]);

                //高斯曲率，平均曲率，主曲率
                KG[i] = a[i] * c[i] - (b[i] * b[i] / 4);
                H[i] = (a[i] + c[i]) / 2;
                k1[i] = H[i] + Math.Sqrt(H[i] * H[i] - KG[i]);
                k2[i] = H[i] - Math.Sqrt(H[i] * H[i] - KG[i]);

                //主方向E1，E2，其中O为计算角度
                O[i] = 0.5 * Math.Asin(b[i] / (k2[i] - k1[i]));
                IVector3D tempE1 = new Vector3DClass();
                IVector3D tranE1 = new Vector3DClass();
                IVector3D momE1 = new Vector3DClass();
                IVector3D briE1 = new Vector3DClass();
                IVector3D tempE2 = new Vector3DClass();
                IVector3D tranE2 = new Vector3DClass();
                IVector3D momE2 = new Vector3DClass();
                IVector3D briE2 = new Vector3DClass();
                if (E1[i].IsEmpty)
                    continue;
                if (E2[i].IsEmpty)
                    continue;
                Material.Clone(E1[i], tempE1);
                Material.Clone(E2[i], tempE2);
                Material.Clone(E1[i], tranE1);
                Material.Clone(E2[i], tranE2);
                Material.Clone(E1[i], momE1);
                Material.Clone(E2[i], momE2);
                Material.Clone(E1[i], briE1);
                Material.Clone(E2[i], briE2);

                //e1
                tempE1.Magnitude = (momE1.Magnitude) * (Math.Cos(O[i]));
                tempE2.Magnitude = (momE2.Magnitude) * (Math.Sin(O[i]));
                e1[i] = new Vector3DClass();
                if (tempE1.IsEmpty)
                    continue;
                if (tempE2.IsEmpty)
                    continue;
                e1[i].ConstructAddVector(tempE1, tempE2);

                //e2
                tranE1.Magnitude = ((briE1.Magnitude) * (Math.Sin(O[i])));
                tranE2.Magnitude = ((briE2.Magnitude) * (Math.Cos(O[i])));
                e2[i] = new Vector3DClass();
                if (tranE1.IsEmpty)
                    continue;
                if (tranE1.IsEmpty)
                    continue;
                e2[i].ConstructSubtractVector(tranE2, tranE1);

                //O数组
                //更新其ID
                int u = PointLocation(Point[i]);
                Oe1[u] = e1[i];
                Oe2[u] = e2[i];
                Ok1[u] = k1[i];
                Ok2[u] = k2[i];
                OG[u] = KG[i];
                OM[u] = H[i];
            }

            //Entropy
            for (int i = 0; i < NodeCount - 4; i++)
            {
                HEntropy[i] = Material.Entropy(ONode[i], OM, tin);
                OME[i] = HEntropy[i];
                GEntropy[i] = Material.Entropy(ONode[i], OG, tin);
                OGE[i] = GEntropy[i];
            }

            //得到每一个节点处的标准差
            //得到揭露Fd的钻孔编号的索引
            int BoreCount;
            double[,] BoreXY;
            string[] XYLine = File.ReadAllLines(exposepath);
            BoreCount = XYLine.Length;
            BoreXY = new double[BoreCount, 2];
            int l = 0;
            for (int j = 0; j < BoreCount; j++)
            {
                string m = XYLine[j];
                string M = new System.Text.RegularExpressions.Regex("[\\s]+").Replace(m, " ");
                string[] mArray = M.Split(' ');
                BoreXY[j, 0] = Convert.ToDouble(mArray[2]);
                BoreXY[j, 1] = Convert.ToDouble(mArray[1]);

                //得到改变A的钻孔位于数组的序号                                                                                                                                       //3个？
                string ChangA = Convert.ToString(mArray[0]);
                //if ((ChangA == CreateNodes.SelectBore1) || (ChangA == CreateNodes.SelectBore2) || (ChangA == CreateNodes.SelectBore3))
                if ((ChangA == CreateNodes.SelectBore1))
                {
                    CreateNodes.boreLocationChangeA[l] = j;
                    l++;
                }
            }

            //未扰动
            if (CreateNodes.pattern == 0)
            {
                for (int i = 0; i < NodeCount; i++)
                {
                    OSTD[i] = CreateNodes.SigmaBase;
                }
            }

           //正太扰动
            else
            {
                //离钻孔的最近距离
                double[] NearestD = new double[NodeCount];
                //保存辐射它的钻孔
                int[] InfluenceBore = new int[NodeCount];
                for (int i = 0; i < NodeCount - 4; i++)
                {
                    NearestD[i] = 10000000;
                    for (int j = 0; j < BoreCount; j++)
                    {
                        double Distance = Math.Sqrt((ONode[i].X - BoreXY[j, 0]) * (ONode[i].X - BoreXY[j, 0]) + (ONode[i].Y - BoreXY[j, 1]) * (ONode[i].Y - BoreXY[j, 1]));
                        if (NearestD[i] > Distance)
                        {
                            NearestD[i] = Distance;
                            //找出辐射它的钻孔
                            InfluenceBore[i] = j;
                        }
                    }
                }

                //得到每个节点的排序编号(由小到大)
                double[] Boreorder = new double[NodeCount];
                for (int i = 0; i < NodeCount - 4; i++)
                {
                    Boreorder[i] = 1;
                    for (int j = 0; j < NodeCount - 4; j++)
                    {
                        if (NearestD[i] > NearestD[j])
                        {
                            Boreorder[i]++;
                        }
                    }
                }

                //排序完后，筛选出三个等间节点序号
                ChooseBore = new int[3];

                for (int k = 0; k < NodeCount - 4; k++)
                {
                    if (Boreorder[k] == 1)
                    {
                        ChooseBore[0] = k;
                    }

                    if (Boreorder[k] == NodeCount - 4)
                    {
                        ChooseBore[2] = k;
                    }
                }

                //每个节点对应的0---1之间的指数模型的y值
                //最大距离差
                double MaxDistance = NearestD[ChooseBore[2]];

                //每个距离减去最小的距离/最大距离减去最小的距离-->代入指数模型得到y值
                for (int i = 0; i < NodeCount; i++)
                {
                    double Ex = NearestD[i] / MaxDistance;

                    if (CreateNodes.pattern == 1)
                    {
                        OSTD[i] = (Efunction(Ex, CreateNodes.NdistPara) + 1) * CreateNodes.SigmaBase;

                        //筛出改A的
                        int nn = CreateNodes.boreLocationChangeA.Length;
                        for (int t = 0; t < nn; t++)
                        {
                            if (InfluenceBore[i] == CreateNodes.boreLocationChangeA[t])
                            {
                                OSTD[i] = (Efunction(Ex, CreateNodes.NdistPara) + 1) * CreateNodes.SigmaChange;
                                break;
                            }
                        }
                    }
                    if (CreateNodes.pattern == 2)
                    {
                        OSTD[i] = (Efunction(Ex, CreateNodes.UdistPara) + 1) * CreateNodes.IntervalBase;
                    }
                }

                //选中间值
                if (CreateNodes.pattern == 1)
                {
                    double median = (CreateNodes.NdistPara + CreateNodes.SigmaBase) / 2;

                    for (int i = 0; i < NodeCount; i++)
                    {
                        if (Math.Abs(OSTD[i] - median) < 0.02)                                                                      //
                        {
                            ChooseBore[1] = i;
                            break;
                        }
                    }
                }

                if (CreateNodes.pattern == 2)
                {
                    double median = (CreateNodes.UdistPara + CreateNodes.IntervalBase) / 2;

                    for (int i = 0; i < NodeCount; i++)
                    {
                        if (Math.Abs(OSTD[i] - median) < 0.02)                                                                      //
                        {
                            ChooseBore[1] = i;
                            break;
                        }
                    }
                }
            }
        }
        #endregion DifferentialGeometry

        #region CreateFeature
        #region Prepare
        //要素类
        public IFeatureClass CreateStandaloneFeatureClass(IWorkspace workspace, String featureClassName,
        IFields fieldsCollection)
        {
            IFeatureWorkspace featureWorkspace = (IFeatureWorkspace)workspace;
            IFeatureClassDescription fcDesc = new FeatureClassDescriptionClass();
            IObjectClassDescription ocDesc = (IObjectClassDescription)fcDesc;

            // Use IFieldChecker to create a validated fields collection.
            IFieldChecker fieldChecker = new FieldCheckerClass();
            IEnumFieldError enumFieldError = null;
            IFields validatedFields = null;
            fieldChecker.ValidateWorkspace = workspace;
            fieldChecker.Validate(fieldsCollection, out enumFieldError, out validatedFields);

            // The enumFieldError enumerator can be inspected at this point to determine 
            // which fields were modified during validation.
            IFeatureClass featureClass = featureWorkspace.CreateFeatureClass(featureClassName, validatedFields,
            ocDesc.InstanceCLSID, ocDesc.ClassExtensionCLSID, esriFeatureType.esriFTSimple, fcDesc.ShapeFieldName, "");
            return featureClass;
        }

        //字段集
        private IFields CreateFieldsCollectionForFeatueClass(ISpatialReference spatialReference,
            esriGeometryType geometryType)
        {
            //Use the feature class description to return the required fields in a fields collection.
            IFeatureClassDescription fcDesc = new FeatureClassDescriptionClass();
            IObjectClassDescription ocDesc = (IObjectClassDescription)fcDesc;

            //Create the fields using the required fields method.
            IFields fields = ocDesc.RequiredFields;

            //Locate the shape field with the name from the feature class description.
            int shapeFieldIndex = fields.FindField(fcDesc.ShapeFieldName);
            IField shapeField = fields.get_Field(shapeFieldIndex);

            //Modify the GeometryDef object before using the fields collection to create a
            //feature class.
            IGeometryDef geometryDef = shapeField.GeometryDef;
            IGeometryDefEdit geometryDefEdit = (IGeometryDefEdit)geometryDef;

            //Alter the feature class geometry type to the type we need.
            geometryDefEdit.GeometryType_2 = geometryType;
            geometryDefEdit.HasZ_2 = true;
            geometryDefEdit.HasM_2 = true;
            geometryDefEdit.GridCount_2 = 1;

            //Set the first grid size to zero and allow ArcGIS to determine a valid grid size.
            geometryDefEdit.set_GridSize(0, 0);
            geometryDefEdit.SpatialReference_2 = spatialReference;

            //If it is a point collection, add Field
            if (geometryType == esriGeometryType.esriGeometryPoint)
            {
                //Create a user-defined double field for "Gauss"
                IFieldsEdit fieldsEdit = (IFieldsEdit)fields;
                IField incomeField = new FieldClass();
                IFieldEdit incomeFieldEdit = (IFieldEdit)incomeField;
                incomeFieldEdit.AliasName_2 = "Gauss";
                incomeFieldEdit.Editable_2 = true;
                incomeFieldEdit.IsNullable_2 = false;
                incomeFieldEdit.Name_2 = "GaussCurvature";
                incomeFieldEdit.Precision_2 = 2;
                incomeFieldEdit.Scale_2 = 5;
                incomeFieldEdit.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit.AddField(incomeField);
                //Create a user-defined double field for "H"
                IFieldsEdit fieldsEdit2 = (IFieldsEdit)fields;
                IField incomeField2 = new FieldClass();
                IFieldEdit incomeFieldEdit2 = (IFieldEdit)incomeField2;
                incomeFieldEdit2.AliasName_2 = "M";
                incomeFieldEdit2.Editable_2 = true;
                incomeFieldEdit2.IsNullable_2 = false;
                incomeFieldEdit2.Name_2 = "Meancurvature";
                incomeFieldEdit2.Precision_2 = 2;
                incomeFieldEdit2.Scale_2 = 5;
                incomeFieldEdit2.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit2.AddField(incomeField2);
                //Create a user-defined double field for "k1"
                IFieldsEdit fieldsEdit3 = (IFieldsEdit)fields;
                IField incomeField3 = new FieldClass();
                IFieldEdit incomeFieldEdit3 = (IFieldEdit)incomeField3;
                incomeFieldEdit3.AliasName_2 = "k1";
                incomeFieldEdit3.Editable_2 = true;
                incomeFieldEdit3.IsNullable_2 = false;
                incomeFieldEdit3.Name_2 = "BigCurvature";
                incomeFieldEdit3.Precision_2 = 2;
                incomeFieldEdit3.Scale_2 = 5;
                incomeFieldEdit3.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit3.AddField(incomeField3);
                //Create a user-defined double field for "k2"
                IFieldsEdit fieldsEdit4 = (IFieldsEdit)fields;
                IField incomeField4 = new FieldClass();
                IFieldEdit incomeFieldEdit4 = (IFieldEdit)incomeField4;
                incomeFieldEdit4.AliasName_2 = "k2";
                incomeFieldEdit4.Editable_2 = true;
                incomeFieldEdit4.IsNullable_2 = false;
                incomeFieldEdit4.Name_2 = "SmallCurvature";
                incomeFieldEdit4.Precision_2 = 2;
                incomeFieldEdit4.Scale_2 = 5;
                incomeFieldEdit4.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit4.AddField(incomeField4);
                //Create a user-defined double field for "X"
                IFieldsEdit fieldsEdit5 = (IFieldsEdit)fields;
                IField incomeField5 = new FieldClass();
                IFieldEdit incomeFieldEdit5 = (IFieldEdit)incomeField5;
                incomeFieldEdit5.AliasName_2 = "X";
                incomeFieldEdit5.Editable_2 = true;
                incomeFieldEdit5.IsNullable_2 = false;
                incomeFieldEdit5.Name_2 = "X";
                incomeFieldEdit5.Precision_2 = 2;
                incomeFieldEdit5.Scale_2 = 5;
                incomeFieldEdit5.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit5.AddField(incomeField5);
                //Create a user-defined double field for "Y"
                IFieldsEdit fieldsEdit6 = (IFieldsEdit)fields;
                IField incomeField6 = new FieldClass();
                IFieldEdit incomeFieldEdit6 = (IFieldEdit)incomeField6;
                incomeFieldEdit6.AliasName_2 = "Y";
                incomeFieldEdit6.Editable_2 = true;
                incomeFieldEdit6.IsNullable_2 = false;
                incomeFieldEdit6.Name_2 = "Y";
                incomeFieldEdit6.Precision_2 = 2;
                incomeFieldEdit6.Scale_2 = 5;
                incomeFieldEdit6.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit6.AddField(incomeField6);
                //Create a user-defined double field for "Z"
                IFieldsEdit fieldsEdit7 = (IFieldsEdit)fields;
                IField incomeField7 = new FieldClass();
                IFieldEdit incomeFieldEdit7 = (IFieldEdit)incomeField7;
                incomeFieldEdit7.AliasName_2 = "Z";
                incomeFieldEdit7.Editable_2 = true;
                incomeFieldEdit7.IsNullable_2 = false;
                incomeFieldEdit7.Name_2 = "Z";
                incomeFieldEdit7.Precision_2 = 2;
                incomeFieldEdit7.Scale_2 = 5;
                incomeFieldEdit7.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit7.AddField(incomeField7);
                //Create a user-defined double field for "ArrayCount"
                IFieldsEdit fieldsEdit8 = (IFieldsEdit)fields;
                IField incomeField8 = new FieldClass();
                IFieldEdit incomeFieldEdit8 = (IFieldEdit)incomeField8;
                incomeFieldEdit8.AliasName_2 = "ArrayCount";
                incomeFieldEdit8.Editable_2 = true;
                incomeFieldEdit8.IsNullable_2 = false;
                incomeFieldEdit8.Name_2 = "ArrayCount";
                incomeFieldEdit8.Precision_2 = 2;
                incomeFieldEdit8.Scale_2 = 5;
                incomeFieldEdit8.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit8.AddField(incomeField8);
                //平均曲率熵
                IFieldsEdit fieldsEdit19 = (IFieldsEdit)fields;
                IField incomeField19 = new FieldClass();
                IFieldEdit incomeFieldEdit19 = (IFieldEdit)incomeField19;
                incomeFieldEdit19.AliasName_2 = "MEntropy";
                incomeFieldEdit19.Editable_2 = true;
                incomeFieldEdit19.IsNullable_2 = false;
                incomeFieldEdit19.Name_2 = "MEntropy";
                incomeFieldEdit19.Precision_2 = 2;
                incomeFieldEdit19.Scale_2 = 5;
                incomeFieldEdit19.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit19.AddField(incomeField19);
                //高斯曲率熵
                IFieldsEdit fieldsEdit119 = (IFieldsEdit)fields;
                IField incomeField119 = new FieldClass();
                IFieldEdit incomeFieldEdit119 = (IFieldEdit)incomeField119;
                incomeFieldEdit119.AliasName_2 = "GEntropy";
                incomeFieldEdit119.Editable_2 = true;
                incomeFieldEdit119.IsNullable_2 = false;
                incomeFieldEdit119.Name_2 = "GEntropy";
                incomeFieldEdit119.Precision_2 = 2;
                incomeFieldEdit119.Scale_2 = 5;
                incomeFieldEdit119.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit119.AddField(incomeField119);
                //标准差
                IFieldsEdit fieldsEdit199 = (IFieldsEdit)fields;
                IField incomeField199 = new FieldClass();
                IFieldEdit incomeFieldEdit199 = (IFieldEdit)incomeField199;
                incomeFieldEdit199.AliasName_2 = "STD";
                incomeFieldEdit199.Editable_2 = true;
                incomeFieldEdit199.IsNullable_2 = false;
                incomeFieldEdit199.Name_2 = "STD";
                incomeFieldEdit199.Precision_2 = 2;
                incomeFieldEdit199.Scale_2 = 5;
                incomeFieldEdit199.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit199.AddField(incomeField199);
            }
            //If it is a polyline, add Field
            if (geometryType == esriGeometryType.esriGeometryPolyline)
            {
                //Create a user-defined double field for "Magnitude"
                IFieldsEdit fieldsEdit = (IFieldsEdit)fields;
                IField incomeField = new FieldClass();
                IFieldEdit incomeFieldEdit = (IFieldEdit)incomeField;
                incomeFieldEdit.AliasName_2 = "VectorMagnitude";
                incomeFieldEdit.Editable_2 = true;
                incomeFieldEdit.IsNullable_2 = false;
                incomeFieldEdit.Name_2 = "VectorMagnitude";
                incomeFieldEdit.Precision_2 = 2;
                incomeFieldEdit.Scale_2 = 5;
                incomeFieldEdit.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit.AddField(incomeField);
            }

            if (geometryType == esriGeometryType.esriGeometryMultiPatch)
            {
                //Create a user-defined double field for "IsRegular"
                IFieldsEdit fieldsEdit3 = (IFieldsEdit)fields;
                IField incomeField3 = new FieldClass();
                IFieldEdit incomeFieldEdit3 = (IFieldEdit)incomeField3;
                incomeFieldEdit3.AliasName_2 = "IsRegular";
                incomeFieldEdit3.Editable_2 = true;
                incomeFieldEdit3.IsNullable_2 = false;
                incomeFieldEdit3.Name_2 = "IsRegular";
                incomeFieldEdit3.Precision_2 = 2;
                incomeFieldEdit3.Scale_2 = 5;
                incomeFieldEdit3.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit3.AddField(incomeField3);
            }
            return fields;
        }

        public void CreateFeatureClass()
        {
            string NodeClass = "Node";
            string NodeVectorLineClassOne = "e1";
            string NodeVectorArrowClassOne = "e1Arrow";
            string NodeVectorLineClassTwo = "e2";
            string NodeVectorArrowClassTwo = "e2Arrow";
            string UnitNormalVectorOfNodeLineClass = "NormalVector";
            string UnitNormalVectorOfNodeArrowClass = "NormalVectorArrow";

            //传递坐标系
            IGeoDataset geoDataset = tin as IGeoDataset;
            ISpatialReference spatialReference = geoDataset.SpatialReference;
            //创建要素类
            NodeFields = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPoint);
            NodeFeatureClass = CreateStandaloneFeatureClass(Workspace, NodeClass, NodeFields);

            UnitNormalVectorOfNodeLinefields = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPolyline);
            UnitNormalVectorOfNodeArrowfields = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryMultiPatch);
            UnitNormalVectorOfNodeLineFeatureClass = CreateStandaloneFeatureClass(Workspace, UnitNormalVectorOfNodeLineClass, UnitNormalVectorOfNodeLinefields);
            UnitNormalVectorOfNodeArrowFeatureClass = CreateStandaloneFeatureClass(Workspace, UnitNormalVectorOfNodeArrowClass, UnitNormalVectorOfNodeArrowfields);

            NodeVectorLinefieldsOne = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPolyline);
            NodeVectorArrowfieldsOne = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryMultiPatch);
            NodeVectorLineFeatureClassOne = CreateStandaloneFeatureClass(Workspace, NodeVectorLineClassOne, NodeVectorLinefieldsOne);
            NodeVectorArrowFeatureClassOne = CreateStandaloneFeatureClass(Workspace, NodeVectorArrowClassOne, NodeVectorArrowfieldsOne);

            NodeVectorLinefieldsTwo = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPolyline);
            NodeVectorArrowfieldsTwo = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryMultiPatch);
            NodeVectorLineFeatureClassTwo = CreateStandaloneFeatureClass(Workspace, NodeVectorLineClassTwo, NodeVectorLinefieldsTwo);
            NodeVectorArrowFeatureClassTwo = CreateStandaloneFeatureClass(Workspace, NodeVectorArrowClassTwo, NodeVectorArrowfieldsTwo);

        }
        #endregion Prepare

        public void CreateFeature()
        {
            //移动放缩矢量箭头
            IVector3D VerticalVector = Material.ConstructVector3D(0, 0, 1);
            IPoint originPoint = Material.ConstructPoint3D(0, 0, 0);
            Material.MakeZAware(originPoint as IGeometry);
            double XScale = VectorModel.XScale;                                                                             //
            double YScale = VectorModel.YScale;                                                                             //    适度调整
            double ZScale = VectorModel.ZScale;                                                                             //

            #region NodeCurvature
            //创建顶点
            for (int i = 4; i < NodeCount; i++)
            {
                if (NormalVectorOfNode[i].IsEmpty)
                    continue;
                IFeature feature1 = NodeFeatureClass.CreateFeature();
                feature1.Shape = Point[i];
                int index1 = feature1.Fields.FindField("GaussCurvature");
                feature1.set_Value(index1, KG[i]);
                int index2 = feature1.Fields.FindField("Meancurvature");
                feature1.set_Value(index2, H[i]);
                int index3 = feature1.Fields.FindField("BigCurvature");
                feature1.set_Value(index3, k1[i]);
                int index4 = feature1.Fields.FindField("SmallCurvature");
                feature1.set_Value(index4, k2[i]);
                int index51 = feature1.Fields.FindField("STD");
                int w = NodeLocation(tin.GetNode(i + 1));
                feature1.set_Value(index51, OSTD[w]);
                int index11 = feature1.Fields.FindField("X");
                feature1.set_Value(index11, Point[i].X);
                int index12 = feature1.Fields.FindField("Y");
                feature1.set_Value(index12, Point[i].Y);
                int index13 = feature1.Fields.FindField("Z");
                feature1.set_Value(index13, Point[i].Z);
                //更新其ID
                int u = PointLocation(Point[i]);
                int index5 = feature1.Fields.FindField("ArrayCount");
                feature1.set_Value(index5, u);
                int index14 = feature1.Fields.FindField("MEntropy");
                feature1.set_Value(index14, HEntropy[u]);
                int index114 = feature1.Fields.FindField("GEntropy");
                feature1.set_Value(index114, GEntropy[u]);
                feature1.Store();
            }
            #endregion NodeCurvature

            #region NodeNormalVector
            //创建基于顶点的法矢量
            for (int i = 4; i < NodeCount; i++)
            {
                if (Point[i].IsEmpty)
                    continue;
                IRay ray = new RayClass();
                ray.Origin = Point[i];
                if (NormalVectorOfNode[i].IsEmpty)
                    continue;
                ray.Vector = NormalVectorOfNode[i];
                IPoint endPoint = new PointClass();
                Material.MakeZAware(endPoint as IGeometry);
                //设置矢量可视化长度
                endPoint = ray.GetPointAtDistance(VectorModel.NLength_D);                                              //  在VM类里面                                                                                                     
                //添加两点构成向量的线端点
                IPointCollection VectorPointCollection = new PolylineClass();
                Material.MakeZAware(VectorPointCollection as IGeometry);
                VectorPointCollection.AddPoint(Point[i], ref _missing, ref _missing);
                VectorPointCollection.AddPoint(endPoint, ref _missing, ref _missing);
                //为向量末端添加箭头符号
                IGeometry Arrow = Material.GetArrow();
                double Inclination = NormalVectorOfNode[i].Inclination;
                //计算与竖向矢量夹角
                double degreesOfRotation = Math.Acos((NormalVectorOfNode[i].DotProduct(VerticalVector)) / ((NormalVectorOfNode[i].Magnitude) * (VerticalVector.Magnitude)));
                ITransform3D transform3D = Arrow as ITransform3D;
                //放缩
                transform3D.Scale3D(originPoint, XScale, YScale, ZScale);
                //转动
                if (degreesOfRotation != 0)
                {
                    double angleOfRotationInRadians = degreesOfRotation;
                    IVector3D axisOfRotationVector3D = new Vector3DClass();
                    axisOfRotationVector3D.ConstructCrossProduct(VerticalVector, NormalVectorOfNode[i]);
                    transform3D.RotateVector3D(axisOfRotationVector3D, angleOfRotationInRadians);
                }
                //平移
                if (endPoint.IsEmpty)
                    continue;
                transform3D.Move3D(endPoint.X - originPoint.X, endPoint.Y - originPoint.Y, endPoint.Z - originPoint.Z);
                //创建要素
                IFeature feature10 = UnitNormalVectorOfNodeLineFeatureClass.CreateFeature();
                feature10.Shape = VectorPointCollection as IGeometry;
                int index1 = feature10.Fields.FindField("VectorMagnitude");
                feature10.set_Value(index1, NormalVectorOfNode[i].Magnitude);
                feature10.Store();
                IFeature feature11 = UnitNormalVectorOfNodeArrowFeatureClass.CreateFeature();
                feature11.Shape = Arrow as IMultiPatch;
                feature11.Store();
            }
            #endregion NodeNormalVector

            #region PrincipalDirection
            //创建基于顶点的主方向矢量线
            //e1
            for (int i = 4; i < NodeCount; i++)
            {
                if (Point[i].IsEmpty)
                    continue;
                IRay ray = new RayClass();
                ray.Origin = Point[i];
                if (e1[i] == null)
                    continue;
                ray.Vector = e1[i];
                IPoint endPoint = new PointClass();
                Material.MakeZAware(endPoint as IGeometry);
                //设置矢量可视化长度
                endPoint = ray.GetPointAtDistance(Math.Abs(k1[i] * VectorModel.ELengh_D));                       //  在VM类里面        
                //添加两点构成向量的线端点
                IPointCollection VectorPointCollection = new PolylineClass();
                Material.MakeZAware(VectorPointCollection as IGeometry);
                VectorPointCollection.AddPoint(Point[i], ref _missing, ref _missing);
                VectorPointCollection.AddPoint(endPoint, ref _missing, ref _missing);
                //为向量末端添加箭头符号
                IGeometry Arrow = Material.GetArrow();
                if (e1[i].IsEmpty)
                    continue;
                double Inclination = e1[i].Inclination;
                double degreesOfRotation = Math.Acos((e1[i].DotProduct(VerticalVector)) / ((e1[i].Magnitude) * (VerticalVector.Magnitude)));
                ITransform3D transform3D = Arrow as ITransform3D;
                //放缩
                transform3D.Scale3D(originPoint, XScale, YScale, ZScale);
                //转动
                if (degreesOfRotation != 0)
                {
                    double angleOfRotationInRadians = degreesOfRotation;
                    IVector3D axisOfRotationVector3D = new Vector3DClass();
                    axisOfRotationVector3D.ConstructCrossProduct(VerticalVector, e1[i]);
                    transform3D.RotateVector3D(axisOfRotationVector3D, angleOfRotationInRadians);
                }
                //平移
                transform3D.Move3D(endPoint.X - originPoint.X, endPoint.Y - originPoint.Y, endPoint.Z - originPoint.Z);
                //创建要素
                IFeature feature6 = NodeVectorLineFeatureClassOne.CreateFeature();
                if (VectorPointCollection == null)
                    continue;
                feature6.Shape = VectorPointCollection as IGeometry;
                int index1 = feature6.Fields.FindField("VectorMagnitude");
                feature6.set_Value(index1, k1[i]);
                feature6.Store();
                IFeature feature7 = NodeVectorArrowFeatureClassOne.CreateFeature();
                feature7.Shape = Arrow as IMultiPatch;
                if (k1[i] >= 0)
                {
                    IsRegular = 1;
                }
                else
                {
                    IsRegular = 0;
                }
                int index01 = feature7.Fields.FindField("IsRegular");
                feature7.set_Value(index01, IsRegular);
                feature7.Store();
            }

            //e2
            for (int i = 4; i < NodeCount; i++)
            {
                if (Point[i].IsEmpty)
                    continue;
                IRay ray = new RayClass();
                ray.Origin = Point[i];
                if (e2[i] == null)
                    continue;
                ray.Vector = e2[i];
                IPoint endPoint = new PointClass();
                Material.MakeZAware(endPoint as IGeometry);
                //设置矢量可视化长度
                endPoint = ray.GetPointAtDistance(Math.Abs(k2[i] * VectorModel.ELengh_D));
                //添加两点构成向量的线端点
                IPointCollection VectorPointCollection = new PolylineClass();
                Material.MakeZAware(VectorPointCollection as IGeometry);
                VectorPointCollection.AddPoint(Point[i], ref _missing, ref _missing);
                VectorPointCollection.AddPoint(endPoint, ref _missing, ref _missing);
                //为向量末端添加箭头符号
                IGeometry Arrow = Material.GetArrow();
                if (e2[i] == null)
                    continue;
                double Inclination = e2[i].Inclination;
                double degreesOfRotation = Math.Acos((e2[i].DotProduct(VerticalVector)) / ((e2[i].Magnitude) * (VerticalVector.Magnitude)));
                ITransform3D transform3D = Arrow as ITransform3D;
                //放缩
                transform3D.Scale3D(originPoint, XScale, YScale, ZScale);
                //转动
                if (degreesOfRotation != 0)
                {
                    double angleOfRotationInRadians = degreesOfRotation;
                    IVector3D axisOfRotationVector3D = new Vector3DClass();
                    axisOfRotationVector3D.ConstructCrossProduct(VerticalVector, e2[i]);
                    transform3D.RotateVector3D(axisOfRotationVector3D, angleOfRotationInRadians);
                }
                //平移
                if (endPoint.IsEmpty)
                    continue;
                transform3D.Move3D(endPoint.X - originPoint.X, endPoint.Y - originPoint.Y, endPoint.Z - originPoint.Z);
                //创建要素
                IFeature feature8 = NodeVectorLineFeatureClassTwo.CreateFeature();
                if (VectorPointCollection == null)
                    continue;
                feature8.Shape = VectorPointCollection as IGeometry;
                int index1 = feature8.Fields.FindField("VectorMagnitude");
                feature8.set_Value(index1, k2[i]);
                feature8.Store();
                IFeature feature9 = NodeVectorArrowFeatureClassTwo.CreateFeature();
                feature9.Shape = Arrow as IMultiPatch;
                if (k2[i] >= 0)
                {
                    IsRegular = 1;
                }
                else
                {
                    IsRegular = 0;
                }
                int index01 = feature9.Fields.FindField("IsRegular");
                feature9.set_Value(index01, IsRegular);
                feature9.Store();
            }
            #endregion PrincipalDirection

            //写入txt文件

            StreamWriter PointsFile = new StreamWriter(pointsFile);
            bool Head = true;
            for (int i = 0; i < NodeCount - 4; i++)
            {
                if (Head == true)
                {
                    PointsFile.WriteLine("GOCAD VSet" + "\r\n"
                        + "HEADER" + "\r\n"
                        + "{name: points}" + "\r\n"
                        + "GOCAD_ORIGINAL_COORDINATE_SYSTEM\nNAME" + "\r\n"
                        + "PROJECTION Unknown" + "\r\n"
                        + "DATUM Unknown" + "\r\n"
                        + "AXIS_NAME X Y Z" + "\r\n"
                        + "AXIS_UNIT m m m" + "\r\n"
                        + "ZPOSITIVE Elevation" + "\r\n"
                        + "END_ORIGINAL_COORDINATE_SYSTEM" + "\r\n"
                        + "PROPERTIES G M ME GE STD" + "\r\n"
                        + "PROP_LEGAL_RANGES **none**  **none** **none**  **none** **none**  **none** **none**  **none** **none**  **none**" + "\r\n"
                        + "NO_DATA_VALUES -99999 -99999 -99999 -99999 -99999" + "\r\n"
                        + "PROPERTY_CLASSES g m me ge std" + "\r\n"
                        + "PROPERTY_KINDS \"Real Number\" \"Real Number\" \"Real Number\" \"Real Number\" \"Real Number\"" + "\r\n"
                        + "PROPERTY_SUBCLASSES QUANTITY Float QUANTITY Float QUANTITY Float QUANTITY Float QUANTITY Float" + "\r\n"
                        + "ESIZES 1  1  1  1  1" + "\r\n"
                        + "UNITS unitless unitless unitless unitless unitless" + "\r\n"
                        + "PROPERTY_CLASS_HEADER X {" + "\r\n" + "kind: X" + "\r\n" + "unit: m" + "\r\n" + "pclip: 99}" + "\r\n"
                        + "PROPERTY_CLASS_HEADER Y {" + "\r\n" + "kind: Y" + "\r\n" + "unit: m" + "\r\n" + "pclip: 99}" + "\r\n"
                        + "PROPERTY_CLASS_HEADER Z {" + "\r\n" + "kind: Depth\nunit: m" + "\r\n" + "is_z: on" + "\r\n" + "pclip: 99}" + "\r\n"
                        + "PROPERTY_CLASS_HEADER g {" + "\r\n" + "kind: Real Number" + "\r\n" + "unit: unitless" + "\r\n" + "pclip: 99}" + "\r\n"
                        + "PROPERTY_CLASS_HEADER m {" + "\r\n" + "kind: Real Number" + "\r\n" + "unit: unitless" + "\r\n" + "pclip: 99}" + "\r\n"
                        + "PROPERTY_CLASS_HEADER me {" + "\r\n" + "kind: Real Number" + "\r\n" + "unit: unitless" + "\r\n" + "pclip: 99}" + "\r\n"
                        + "PROPERTY_CLASS_HEADER ge {" + "\r\n" + "kind: Real Number" + "\r\n" + "unit: unitless" + "\r\n" + "pclip: 99}" + "\r\n"
                        + "PROPERTY_CLASS_HEADER std {" + "\r\n" + "kind: Real Number" + "\r\n" + "unit: unitless" + "\r\n" + "pclip: 99}" + "\r\n"
                        );
                    Head = false;
                }

                string vsX = (OPoint[i].X).ToString();
                string vsY = (OPoint[i].Y).ToString();
                string vsZ = (OPoint[i].Z).ToString();
                string vsG = (OG[i]).ToString();
                string vsM = (OM[i]).ToString();
                string vsME = (OME[i]).ToString();
                string vsGE = (OGE[i]).ToString();
                string vsSTD = (OSTD[i]).ToString();

                PointsFile.Write("PVRTX " + i + " ");
                PointsFile.Write(vsX + " ");
                PointsFile.Write(vsY + " ");
                PointsFile.Write(vsZ + " ");
                PointsFile.Write(vsG + " ");
                PointsFile.Write(vsM + " ");
                PointsFile.Write(vsME + " ");
                PointsFile.Write(vsGE + " ");
                PointsFile.WriteLine(vsSTD + " ");
            }
            PointsFile.Write("END");
            PointsFile.Close();
        }
        #region 排序函数
        public int NodeLocation(ITinNode NodeL)
        {
            int v = (int)((NodeL.X - (XMin - (RasterSize / 2))) / RasterSize) + RowCount * (int)((NodeL.Y - (YMin - (RasterSize / 2))) / (RasterSize));
            return v;
        }

        public int PointLocation(IPoint PointL)
        {
            int v = (int)((PointL.X - (XMin - (RasterSize / 2))) / RasterSize) + RowCount * (int)((PointL.Y - (YMin - (RasterSize / 2))) / (RasterSize));
            return v;
        }

        //指数模型
        public double Efunction(double x, double range)
        {
            double Y = (Math.Exp(x) - 1) / (Math.Exp(1) - 1) * range;
            return Y;
        }
        #endregion
        #endregion CreateFeature
    }
}



