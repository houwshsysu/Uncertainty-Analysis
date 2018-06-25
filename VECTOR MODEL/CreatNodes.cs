using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using ESRI.ArcGIS;
using ESRI.ArcGIS.Geodatabase;
using ESRI.ArcGIS.Geometry;
using System.Collections;

/* >>>Create the nodes of each isosurface
 * 
 * 
 * 
 * 
 * 
 * 
 */

namespace VECTOR_MODEL
{
    class CreateNodes
    {
        //---------------------------------------------------------------------------------------------
        //how to disturbed
        //正态分布-->1        均匀分布-->2
        public static double pattern = 1;
        //public static double pattern = 2;

        //扰动变量
        public static double NdistPara = 3;                              // A
        //正太分布误差分布函数参数                           
        public static double SigmaBase = 1;                            // sigma
        public static double SigmaChange = 1.5;                     // one change sigma

        //均匀分布误差分布函数参数
        public static double IntervalBase = 0.91;
        public static double UdistPara = 2;

        string workspacepath = @"C:\Users\B410\Desktop\beforeSchool\optimization\FdISOsurfaceNode3-1-1.5.gdb";
        string tinpath = @"C:\Users\B410\Desktop\C&G\Data\fd";
        string boreholepath = @"C:\Users\B410\Desktop\beforeSchool\BBBB\Borehole3-1-1.5.txt";
        //揭露处的钻孔XY 
        string XYpath = @"C:\Users\B410\Desktop\C&G\Data\briefBore\expose.txt";
        string OutTinPath = @"C:\Users\B410\Desktop\beforeSchool\optimization\FdISOsurface3-1-1.5\tin_node";

        //修改A的钻孔：01
        public static string SelectBore1 = "MKZ3-YJL-01";
        public static int[] boreLocationChangeA = new int[1];

        //每一列点的稠密程度
        int density = 200;                                               //1m*个点
        double HalfBoreLength = 20;                             //上下各*m （长度不够会抛出异常）
        int SurfaceCount = 50;                                       //需要的曲面数
        double RasterSize = 5;                                       //栅格大小
        int SimulationCounts = 10000;                          //模拟次数
        //---------------------------------------------------------------------------------------------

        //每个节点的误差分布函数参数
        double[] ENparameters;
        int NodeCount;
        ITinAdvanced2 tin;
        int RowCount;
        double XMin;
        double XMax;
        double YMin;
        double YMax;
        IPoint[] OSPoints;
        ITinNode[] OSNodes;
        IWorkspace Workspace;
        //单列密集点
        IPoint[] DensePoints;
        //50个曲面节点集合
        IPoint[,] SinglePoints;
        //所有节点的密集P
        double[,] DenseProbability;
        //单列曲面节点及其地质属性概率
        double[,] SingleProbability;
        double[,] RegardP;
        IPoint[,] SurfacePoints;
        //三个虚拟的P柱
        IPoint[,] Pcolumn;
        double[,] SurfaceProbability;
        double[,] SurfaceRegardP;
        int DenseCount;
        double[,] SimulationZ;
        StreamWriter Borehole;
        //每个节点的分布函数参数值
        int BoreCount;
        int[] ChooseBoreForDistance;
        int[] ChooseBoreForEN;

        //揭露Fd的钻孔编号的索引
        double[,] BoreXY;

        private static object _missing = Type.Missing;

        public void CreateNode()
        {
            Type factoryType = Type.GetTypeFromProgID("esriDataSourcesGDB.FileGDBWorkspaceFactory");
            IWorkspaceFactory workspaceFactory = (IWorkspaceFactory)Activator.CreateInstance(factoryType);
            Workspace = workspaceFactory.OpenFromFile(workspacepath, 0);
            // Open an existing TIN, edit the following path as appropriate.
            tin = new TinClass();
            tin.Init(tinpath);
            Borehole = new StreamWriter(boreholepath);
            //得到和tin网直接相关的数据
            IEnvelope Envelope = tin.Extent;
            XMin = Envelope.XMin;
            XMax = Envelope.XMax;
            YMin = Envelope.YMin;
            YMax = Envelope.YMax;
            //列数
            RowCount = (int)((XMax - XMin + 0.5 * RasterSize) / RasterSize) + 1;
            NodeCount = tin.NodeCount;
            //钻孔总点数
            DenseCount = (int)HalfBoreLength * 2 * density;

            //点间距
            Double Dstep = 1.0 / density;
            OSNodes = new ITinNode[NodeCount];
            OSPoints = new IPoint[NodeCount];
            for (int i = 4; i < NodeCount; i++)
            {
                //得到按序排列的初始曲面节点  
                ITinNode OSNode = tin.GetNode(i + 1);
                int v = (int)((OSNode.X - (XMin - (RasterSize / 2))) / RasterSize) + RowCount * (int)((OSNode.Y - (YMin - (RasterSize / 2))) / (RasterSize));
                OSNodes[v] = OSNode;
                OSPoints[v] = new PointClass();
                MakeZAware(OSPoints[v]);
                OSNode.QueryAsPoint(OSPoints[v]);
            }

            //一根中心在原点的密集点柱
            double Zinitial = 0;
            double Zlow = Zinitial - HalfBoreLength;
            DensePoints = new IPoint[DenseCount];
            for (int j = 0; j < DenseCount; j++)
            {
                double thisZ = Zlow + j * Dstep;
                IPoint OneDensePoint = new PointClass();
                MakeZAware(OneDensePoint);
                OneDensePoint.X = 0;
                OneDensePoint.Y = 0;
                OneDensePoint.Z = thisZ;
                DensePoints[j] = OneDensePoint;
            }

            //得到揭露Fd的钻孔编号的索引
            string[] XYLine = File.ReadAllLines(XYpath);
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
                string ChangA = Convert.ToString(mArray[0]);
                //if ((ChangA == SelectBore1) || (ChangA == SelectBore2) || (ChangA == SelectBore3))
                if ((ChangA == SelectBore1))
                {
                    boreLocationChangeA[l] = j;
                    l++;
                }
            }
        }

        public void GetEveryKindSingleProbability()
        {
            SimulationZ = new double[NodeCount, SimulationCounts];
            DenseProbability = new double[NodeCount, DenseCount];
            ENparameters = new double[NodeCount];

            //正太扰动
            //离钻孔的最近距离
            double[] NearestD = new double[NodeCount];
            //保存辐射它的钻孔
            int[] InfluenceBore = new int[NodeCount];
            for (int i = 0; i < NodeCount - 4; i++)
            {
                NearestD[i] = 10000000;
                for (int j = 0; j < BoreCount; j++)
                {
                    double Distance = Math.Sqrt((OSNodes[i].X - BoreXY[j, 0]) * (OSNodes[i].X - BoreXY[j, 0]) + (OSNodes[i].Y - BoreXY[j, 1]) * (OSNodes[i].Y - BoreXY[j, 1]));
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
            ChooseBoreForDistance = new int[3];
            for (int k = 0; k < NodeCount - 4; k++)
            {
                if (Boreorder[k] == 1)
                {
                    ChooseBoreForDistance[0] = k;
                }

                if (Boreorder[k] == NodeCount - 4)
                {
                    ChooseBoreForDistance[2] = k;
                }
            }

            //每个节点对应的0---1之间的指数模型的y值
            //最大距离差
            double MaxDistance = NearestD[ChooseBoreForDistance[2]];

            //每个距离/最大距离-->代入指数模型得到y值
            for (int i = 0; i < NodeCount; i++)
            {
                double Ex = NearestD[i] / MaxDistance;

                if (pattern == 1)
                {
                    ENparameters[i] = (Efunction(Ex, NdistPara) + 1) * SigmaBase;
                    //筛出改Sigma的
                    int nn = boreLocationChangeA.Length;
                    for (int t = 0; t < nn; t++)
                    {
                        if (InfluenceBore[i] == boreLocationChangeA[t])
                        {
                            ENparameters[i] = (Efunction(Ex, NdistPara) + 1) * SigmaChange;
                            break;
                        }
                    }
                }

                if (pattern == 2)
                {
                    ENparameters[i] = (Efunction(Ex, UdistPara) + 1) * IntervalBase;
                }
            }

            //根据EN值得到每个节点的排序编号(由小到大)
            double[] BoreorderBaseEN = new double[NodeCount];
            for (int i = 0; i < NodeCount - 4; i++)
            {
                BoreorderBaseEN[i] = 1;
                for (int j = 0; j < NodeCount - 4; j++)
                {
                    if (ENparameters[i] > ENparameters[j])
                    {
                        BoreorderBaseEN[i]++;
                    }
                }
            }

            //排序完后，筛选出三个等间节点序号
            ChooseBoreForEN = new int[3];
            int M = (NodeCount - 4) / 2;
            for (int k = 0; k < NodeCount - 4; k++)
            {
                if (BoreorderBaseEN[k] == 1)
                {
                    ChooseBoreForEN[0] = k;
                }

                if (BoreorderBaseEN[k] == M)
                {
                    ChooseBoreForEN[0] = k;
                }

                if (BoreorderBaseEN[k] == NodeCount - 4)
                {
                    ChooseBoreForEN[2] = k;
                }
            }

            //扰动情形下的
            if (pattern == 1)
            {
                //得到标准正太分布随机数
                double[] NormalStZ = StNormal(SimulationCounts);
                //所有点标准正态分布映射为非标准
                for (int i = 0; i < NodeCount; i++)
                {
                    for (int j = 0; j < SimulationCounts; j++)
                    {
                        SimulationZ[i, j] = NormalStZ[j] * ENparameters[i];
                    }
                }
            }

            if (pattern == 2)
            {
                for (int i = 0; i < NodeCount; i++)
                {
                    double[] UZ = Uniform(SimulationCounts, ENparameters[i]);
                    for (int j = 0; j < SimulationCounts; j++)
                    {
                        SimulationZ[i, j] = UZ[j];
                    }
                }
            }

            //每一节点处的扰动密集点P
            for (int i = 0; i < NodeCount; i++)
            {
                for (int j = 0; j < DenseCount; j++)
                {
                    if ((DensePoints[j]) == null)
                        continue;
                    double IsBed = 0;
                    for (int k = 0; k < SimulationCounts; k++)
                    {
                        if (DensePoints[j].Z < SimulationZ[i, k])
                        {
                            IsBed++;
                        }
                    }
                    DenseProbability[i, j] = (double)IsBed / SimulationCounts;
                }
            }

            //得到需要的单列点及其地质属性概率
            SinglePoints = new IPoint[NodeCount, SurfaceCount];
            SingleProbability = new double[NodeCount, SurfaceCount];
            RegardP = new double[NodeCount, SurfaceCount];
            //0.05---0.95
            double Pstep = 0.9 / (SurfaceCount - 1);
            for (int i = 0; i < NodeCount; i++)
            {
                int Start = -1;
                double Pvalue = 0.95;
                for (int j = 0; j < SurfaceCount; j++)
                {
                    for (int k = Start + 1; k < DenseCount; k++)
                    {
                        if (DenseProbability[i, k] <= Pvalue)
                        {
                            //由下往上
                            SinglePoints[i, j] = DensePoints[k];
                            SingleProbability[i, j] = DenseProbability[i, k];
                            RegardP[i, j] = Pvalue;
                            Start = k;
                            Pvalue = Pvalue - Pstep;
                            break;
                        }
                    }
                }
            }
        }

        //平移得到所有的等值面点及其地质属性概率
        public void GetSurfaceNode()
        {
            IPoint originPoint = Material.ConstructPoint3D(0, 0, 0);
            SurfacePoints = new IPoint[NodeCount, SurfaceCount];
            Pcolumn = new IPoint[3, DenseCount];
            SurfaceProbability = new double[NodeCount, SurfaceCount];
            SurfaceRegardP = new double[NodeCount, SurfaceCount];
            int CV = 0;
            for (int i = 0; i < NodeCount - 4; i++)
            {
                for (int j = 0; j < SurfaceCount; j++)
                {
                    //复制一遍好平移
                    SurfacePoints[i, j] = new PointClass();
                    MakeZAware(SurfacePoints[i, j]);
                    SurfacePoints[i, j].X = SinglePoints[i, j].X;
                    SurfacePoints[i, j].Y = SinglePoints[i, j].Y;
                    SurfacePoints[i, j].Z = SinglePoints[i, j].Z;
                    SurfaceProbability[i, j] = SingleProbability[i, j];
                    SurfaceRegardP[i, j] = RegardP[i, j];
                    //平移
                    ITransform3D Transform3D = SurfacePoints[i, j] as ITransform3D;
                    Transform3D.Move3D(OSPoints[i].X - originPoint.X, OSPoints[i].Y - originPoint.Y, OSPoints[i].Z - originPoint.Z);
                }

                //平移色阶钻孔(编号)
                if (pattern != 0)
                {
                    if ((i == ChooseBoreForEN[0]) || (i == ChooseBoreForEN[1]) || (i == ChooseBoreForEN[2]))
                    {
                        for (int j = 0; j < DenseCount; j++)
                        {
                            Pcolumn[CV, j] = new PointClass();
                            MakeZAware(Pcolumn[CV, j]);
                            Pcolumn[CV, j].X = DensePoints[j].X;
                            Pcolumn[CV, j].Y = DensePoints[j].Y;
                            Pcolumn[CV, j].Z = DensePoints[j].Z;
                            //平移
                            ITransform3D Transform3D = Pcolumn[CV, j] as ITransform3D;
                            Transform3D.Move3D(OSPoints[i].X - originPoint.X, OSPoints[i].Y - originPoint.Y, OSPoints[i].Z - originPoint.Z);
                        }
                        CV++;
                    }
                }
            }
        }

        //创建点要素类
        public void CreatFeatureClasses()
        {
            //传递坐标系
            IGeoDataset geoDataset = tin as IGeoDataset;
            ISpatialReference spatialReference = geoDataset.SpatialReference;
            //创建地层节点要素
            for (int i = 0; i < SurfaceCount+1; i++)
            {
                string NodeClass = "Node_" + i;
                IFields Nodefields = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPoint, NodeClass);
                IFeatureClass NodeFeatureClass = CreateStandaloneFeatureClass(Workspace, NodeClass, Nodefields);
                //除掉四个空点
                for (int j = 0; j < NodeCount - 4; j++)
                {
                   IFeature feature = NodeFeatureClass.CreateFeature();
                    //加入初始曲面
                    if (i == SurfaceCount)
                    {
                        feature.Shape = OSPoints[j];
                        int index = feature.Fields.FindField("probability");
                        feature.set_Value(index, 0.5);
                        int index1 = feature.Fields.FindField("RegardP");
                        feature.set_Value(index1, 0.5);
                        feature.Store();
                    }
                    else
                    {
                        feature.Shape = SurfacePoints[j, i];
                        int index = feature.Fields.FindField("probability");
                        feature.set_Value(index, SurfaceProbability[j, i]);
                        int index1 = feature.Fields.FindField("RegardP");
                        feature.set_Value(index1, SurfaceRegardP[j, i]);
                        feature.Store();
                    }
                }
            }

            //创建钻孔点要素
            string NodeCloudClass = "NodeCloud";
            IFields NodeCloudfields = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPoint, NodeCloudClass);
            IFeatureClass NodeCloudFeatureClass = CreateStandaloneFeatureClass(Workspace, NodeCloudClass, NodeCloudfields);

            int CCVV = 0;
            //写钻孔文件
            bool Head = true;

            for (int i = 0; i < NodeCount; i++)
            {
                if ((i == ChooseBoreForEN[0]) || (i == ChooseBoreForEN[1]) || (i == ChooseBoreForEN[2]))
                {
                    //选择钻孔编号
                    for (int j = 0; j < DenseCount; j++)
                    {
                        IFeature feature = NodeCloudFeatureClass.CreateFeature();
                        feature.Shape = Pcolumn[CCVV, j];
                        int index = feature.Fields.FindField("probability");
                        feature.set_Value(index, DenseProbability[i, j]);
                        feature.Store();
                    }

                    for (int v = 0; v < DenseCount; v++)
                    {
                        if (Head == true)
                        {
                            Borehole.WriteLine("GOCAD VSet" + "\r\n"
                                + "HEADER" + "\r\n"
                                + "{name: Borehole}" + "\r\n"
                                + "GOCAD_ORIGINAL_COORDINATE_SYSTEM\nNAME" + "\r\n"
                                + "PROJECTION Unknown" + "\r\n"
                                + "DATUM Unknown" + "\r\n"
                                + "AXIS_NAME X Y Z" + "\r\n"
                                + "AXIS_UNIT m m m" + "\r\n"
                                + "ZPOSITIVE Elevation" + "\r\n"
                                + "END_ORIGINAL_COORDINATE_SYSTEM" + "\r\n"
                                + "PROPERTIES P" + "\r\n"
                                + "PROP_LEGAL_RANGES **none**  **none**" + "\r\n"
                                + "NO_DATA_VALUES -99999" + "\r\n"
                                + "PROPERTY_CLASSES p" + "\r\n"
                                + "PROPERTY_KINDS \"Real Number\"" + "\r\n"
                                + "PROPERTY_SUBCLASSES QUANTITY Float" + "\r\n"
                                + "ESIZES 1" + "\r\n"
                                + "UNITS unitless" + "\r\n"
                                + "PROPERTY_CLASS_HEADER X {" + "\r\n" + "kind: X" + "\r\n" + "unit: m" + "\r\n" + "pclip: 99}" + "\r\n"
                                + "PROPERTY_CLASS_HEADER Y {" + "\r\n" + "kind: Y" + "\r\n" + "unit: m" + "\r\n" + "pclip: 99}" + "\r\n"
                                + "PROPERTY_CLASS_HEADER Z {" + "\r\n" + "kind: Depth\nunit: m" + "\r\n" + "is_z: on" + "\r\n" + "pclip: 99}" + "\r\n"
                                + "PROPERTY_CLASS_HEADER p {" + "\r\n" + "kind: Real Number" + "\r\n" + "unit: unitless" + "\r\n" + "pclip: 99}" + "\r\n"
                                );
                            Head = false;
                        }
                        Borehole.Write("PVRTX " + v + " ");
                        Borehole.Write(Pcolumn[CCVV, v].X + " ");
                        Borehole.Write(Pcolumn[CCVV, v].Y + " ");
                        Borehole.Write(Pcolumn[CCVV, v].Z + " ");

                        //正态分布色阶棒
                        Borehole.WriteLine(DenseProbability[i, v] + " ");

                        //均匀分布色阶棒
                        //if (DenseProbability[i, v] > 0.05 && DenseProbability[i, v] < 0.95)
                        //{
                        //    Borehole.WriteLine(1 + " ");
                        //}
                        //else
                        //{
                        //    Borehole.WriteLine(0 + " ");
                        //}              
                    }
                    CCVV++;
                }
            }
            Borehole.Write("END");
            Borehole.Close();
        }


        //根据创建的点要素做tin
        public void CreatTINs()
        {
            Type factoryType = Type.GetTypeFromProgID("esriDataSourcesGDB.FileGDBWorkspaceFactory");
            IWorkspaceFactory workspaceFactory = (IWorkspaceFactory)Activator.CreateInstance(factoryType);
            Workspace = workspaceFactory.OpenFromFile(workspacepath, 0);
            IFeatureWorkspace featureWorkspace = (IFeatureWorkspace)Workspace;
       
            //new tin 的范围
            tin = new TinClass();
            tin.Init(tinpath);
            IEnvelope EnvO = tin.Extent;
            IEnvelope Env = new EnvelopeClass();
            Env.XMax = EnvO.XMax + 10;
            Env.YMax = EnvO.YMax + 10;
            Env.ZMax = EnvO.ZMax + 100;
            Env.XMin = EnvO.XMin + 10;
            Env.YMin = EnvO.YMin + 10;
            Env.ZMin = EnvO.ZMin + 100;
     
            // Instantiate a new empty TIN.
            ITinEdit[] TinEdit = new ITinEdit[SurfaceCount + 1];
            object overwrite = true;
            for (int i = 0; i < SurfaceCount + 1; i++)
            {
                TinEdit[i] = new TinClass();
                TinEdit[i].InitNew(Env);
                TinEdit[i].SaveAs(OutTinPath + "_" + i, ref overwrite);
            }

            IFeatureClass[] ISOpointFeatureClass = new IFeatureClass[SurfaceCount + 1];
            for (int i = 0; i < SurfaceCount + 1; i++)
            {
                ISOpointFeatureClass[i] = featureWorkspace.OpenFeatureClass("Node_" + i);
                IGeometryCollection MultipointGeometryCollection = new MultipointClass();
                MakeZAware(MultipointGeometryCollection as IGeometry);
                for (int p = 0; p < NodeCount-4; p++)
                {
                    IPoint onePoint = ISOpointFeatureClass[i].GetFeature(p+1).Shape as IPoint;
                    MakeZAware(onePoint);
                    MultipointGeometryCollection.AddGeometry(onePoint);
                }
                (TinEdit[i] as ITinEdit2).SetToConstrainedDelaunay();
                TinEdit[i].AddShapeZ(MultipointGeometryCollection as IGeometry, esriTinSurfaceType.esriTinMassPoint, 0, _missing);
                TinEdit[i].Save();
            }
        }

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
            esriGeometryType geometryType, string className)
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
            //添加地质属性概率字段
            if (geometryType == esriGeometryType.esriGeometryPoint)
            {
                //Create a user-defined double field for "probability"
                IFieldsEdit fieldsEdit = (IFieldsEdit)fields;
                IField incomeField = new FieldClass();
                IFieldEdit incomeFieldEdit = (IFieldEdit)incomeField;
                incomeFieldEdit.AliasName_2 = "probability";
                incomeFieldEdit.Editable_2 = true;
                incomeFieldEdit.IsNullable_2 = false;
                incomeFieldEdit.Name_2 = "probability";
                incomeFieldEdit.Precision_2 = 2;
                incomeFieldEdit.Scale_2 = 5;
                incomeFieldEdit.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit.AddField(incomeField);
                //Create a user-defined double field for "RegardP"
                IFieldsEdit fieldsEdit1 = (IFieldsEdit)fields;
                IField incomeField1 = new FieldClass();
                IFieldEdit incomeFieldEdit1 = (IFieldEdit)incomeField1;
                incomeFieldEdit1.AliasName_2 = "RegardP";
                incomeFieldEdit1.Editable_2 = true;
                incomeFieldEdit1.IsNullable_2 = false;
                incomeFieldEdit1.Name_2 = "RegardP";
                incomeFieldEdit1.Precision_2 = 2;
                incomeFieldEdit1.Scale_2 = 5;
                incomeFieldEdit1.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit1.AddField(incomeField1);
            }
            return fields;
        }
        public static void MakeZAware(IGeometry geometry)
        {
            IZAware zAware = geometry as IZAware;
            zAware.ZAware = true;
        }

        //标准正态分布Box_Muller
        public double[] StNormal(int SimulationTimes)
        {
            Random rand = new Random();
            double[] StZ = new double[SimulationTimes];
            double[] U1 = new double[SimulationTimes];
            double[] U2 = new double[SimulationTimes];
            for (int i = 0; i < SimulationTimes; i++)
            {
                U1[i] = rand.NextDouble();
                U2[i] = rand.NextDouble();
                StZ[i] = Math.Sqrt(-2 * Math.Log(U1[i])) * Math.Cos(2 * Math.PI * U2[i]);
            }
            //返回模拟得到的，用以比较的随机数组
            return StZ;
        }

        //均匀分布(range=半个置信区间的长度)
        private double[] Uniform(int SimulationCount, double HalfRange)
        {
            Random rand = new Random();
            double[] U = new double[SimulationCount];
            double[] UZ = new double[SimulationCount];
            for (int i = 0; i < SimulationCount; i++)
            {
                U[i] = rand.NextDouble();
                UZ[i] = 2 * U[i] * HalfRange - HalfRange;
            }
            return UZ;
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
    }
}