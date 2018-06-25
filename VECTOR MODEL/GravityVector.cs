
/* 得到得到等值面包络体的体积，重心，和重力矢量
 */

using System;
using System.Collections;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ESRI.ArcGIS.Geodatabase;
using ESRI.ArcGIS.Geometry;
using System.IO;

namespace VECTOR_MODEL
{
    class GravityVector
    {
        //---------------------------------------------------------------------------------------------
        //延伸出去的长度
        double AddBoundry = 20;
        //一列的点数（控制精度）
        double ColumnCount = 20000;
        //网格尺寸
        double RasterSize = 5;

        //结果输出的工作空间
        string workSpacePath = @"C:\Users\B410\Desktop\optimization\centroid3.gdb";
        //上底面路径集合
        string UpSurfacePath = @"C:\Users\B410\Desktop\optimization\FuISOsurface1-3\tin_node_";
        //下底面路径集合
        string DownSurfacePath = @"C:\Users\B410\Desktop\optimization\FdISOsurface1-3\tin_node_";
        //最高上底面路径
        string maxUpSurfacePath = @"C:\Users\B410\Desktop\optimization\FuISOsurface1-3\tin_node_49";
        //最低下底面路径
        string minDownSurfacePath = @"C:\Users\B410\Desktop\optimization\FdISOsurface1-3\tin_node_0";

        string Cpath = @"C:\Users\B410\Desktop\VC\Centroid3.";

        string Vpath = @"C:\Users\B410\Desktop\VC\Volume3.txt";

        //箭头放缩
        double XScale = 0.0002;
        double YScale = 0.0002;
        double ZScale = 0.0002;

        //矢量线放缩
        double gLength = 0.0002;
        //---------------------------------------------------------------------------------------------

        IWorkspace Workspace;
        ITinAdvanced2 tinU;
        ITinAdvanced2 tinD;
        ITinAdvanced2 tinUmax;
        ITinAdvanced2 tinDmin;

        double XMin;
        double XMax;
        double YMin;
        double YMax;
        double ZMax;
        double ZMin;
        int RowCount;
        double Llength;
        double Wlength;
        double Hlength;
        double BoxVolume;
        int NodeCount;

        ITinNode[] orderNodeU;
        ITinNode[] orderNodeD;

        StreamWriter Ctxt;
        StreamWriter Vtxt;

        double[] VOLUME = new double[VectorModel.TinCount];
        IPoint[] BARYCENTRE = new IPoint[VectorModel.TinCount];
        IGeometryCollection[] multiPatchGeometryCollection = new IGeometryCollection[VectorModel.TinCount];
        int TINID;
        private static object _missing = Type.Missing;

        //要素类声明
        IFields BarycentreFields;
        IFeatureClass BarycentreFeatureClass;
        IFields GravityVectorFields;
        IFeatureClass GravityVectorFeatureClass;
        IFields ArrowFields;
        IFeatureClass ArrowFeatureClass;
        IFields EncloseFields;
        IFeatureClass[] EncloseFeatureClass = new IFeatureClass[VectorModel.TinCount];

        public void GetTinData(int TinID)
        {
            Type factoryType = Type.GetTypeFromProgID("esriDataSourcesGDB.FileGDBWorkspaceFactory");
            IWorkspaceFactory workspaceFactory = (IWorkspaceFactory)Activator.CreateInstance(factoryType);
            Workspace = workspaceFactory.OpenFromFile(workSpacePath, 0);
            // Open an existing TIN, edit the following path as appropriate.
            tinU = new TinClass();
            tinU.Init(UpSurfacePath + TinID);
            tinD = new TinClass();
            tinD.Init(DownSurfacePath + TinID);
            tinUmax = new TinClass();
            tinUmax.Init(maxUpSurfacePath);
            tinDmin = new TinClass();
            tinDmin.Init(minDownSurfacePath);
            IEnvelope EnvelopeU = tinUmax.Extent;
            IEnvelope EnvelopeD = tinDmin.Extent;
            double MZMax = EnvelopeU.ZMax;
            double MZMin = EnvelopeD.ZMin;
            XMin = EnvelopeU.XMin;
            XMax = EnvelopeU.XMax;
            YMin = EnvelopeU.YMin;
            YMax = EnvelopeU.YMax;
            ZMax = MZMax + AddBoundry;
            ZMin = MZMin - AddBoundry;
            RowCount = (int)((XMax - XMin + 0.5 * RasterSize) / RasterSize) + 1;
            Llength = XMax - XMin;
            Wlength = YMax - YMin;
            Hlength = ZMax - ZMin;
            BoxVolume = Llength * Wlength * Hlength;
            NodeCount = tinDmin.NodeCount;
            orderNodeU = new ITinNode[NodeCount];
            orderNodeD = new ITinNode[NodeCount];
            TINID = TinID;
        }

        public void GetVolume()
        {
            for (int i = 4; i < NodeCount; i++)
            {
                ITinNode NodeU = tinU.GetNode(i + 1);
                int locationU = NodeLocation(NodeU);
                orderNodeU[locationU] = NodeU;
                ITinNode NodeD = tinD.GetNode(i + 1);
                int locationD = NodeLocation(NodeD);
                orderNodeD[locationD] = NodeD;
            }

            double[] SingleSize = new double[NodeCount - 4];
            for (int i = 0; i < NodeCount - 4; i++)
            {
                double Dvalue = orderNodeU[i].Z - orderNodeD[i].Z;
                SingleSize[i] = Dvalue;
            }
            double AllSize = 0;
            for (int i = 0; i < NodeCount - 4; i++)
            {
                AllSize = AllSize + SingleSize[i];
            }
            double BoxSize = (double)(NodeCount - 4) * (Hlength);
            double ratio = AllSize / BoxSize;
            VOLUME[TINID] = BoxVolume * ratio;
        }

        public void GetBarycentre()
        {
            double PointLag = Hlength / (double)(ColumnCount - 1.0);
            double[] upZ = new double[NodeCount - 4];
            double[] downZ = new double[NodeCount - 4];
            double[] DiscreteP = new double[(int)ColumnCount];
            for (int i = 0; i < NodeCount - 4; i++)
            {
                upZ[i] = orderNodeU[i].Z;
                downZ[i] = orderNodeD[i].Z;
            }
            for (int j = 0; j < ColumnCount; j++)
            {
                DiscreteP[j] = ZMin + j * PointLag;
            }
            //比较
            double TotalPCount = 0;
            double TotalX = 0;
            double TotalY = 0;
            double TotalZ = 0;
            for (int i = 0; i < NodeCount - 4; i++)
            {
                for (int j = 0; j < ColumnCount; j++)
                {
                    if ((DiscreteP[j] > downZ[i]) && (DiscreteP[j] < upZ[i]))
                    {
                        TotalX = TotalX + orderNodeD[i].X;
                        TotalY = TotalY + orderNodeD[i].Y;
                        TotalZ = TotalZ + DiscreteP[j];
                        TotalPCount++;
                    }
                }
            }
            //得到重心坐标
            BARYCENTRE[TINID] = new PointClass();
            Material.MakeZAware(BARYCENTRE[TINID]);
            BARYCENTRE[TINID].X = TotalX / TotalPCount;
            BARYCENTRE[TINID].Y = TotalY / TotalPCount;
            BARYCENTRE[TINID].Z = TotalZ / TotalPCount;
        }

        public void CreatEnclose()
        {
            //边界索引
            //UP
            ArrayList[] UpIndexList = new ArrayList[4];
            //DOWN
            ArrayList[] DownIndexList = new ArrayList[4];
            for (int i = 0; i < 4; i++)
            {
                UpIndexList[i] = new ArrayList();
                DownIndexList[i] = new ArrayList();
            }

            //四壁
            for (int i = 0; i < NodeCount - 4; i++)
            {
                if (orderNodeU[i].Y < YMin + RasterSize / 2)
                {
                    UpIndexList[0].Add(i);
                }
                if (orderNodeU[i].X > XMax - RasterSize / 2)
                {
                    UpIndexList[1].Add(i);
                }
                if (orderNodeU[i].Y > YMax - RasterSize / 2)
                {
                    UpIndexList[2].Add(i);
                }
                if (orderNodeU[i].X < XMin + RasterSize / 2)
                {
                    UpIndexList[3].Add(i);
                }
            }

            for (int i = 0; i < 4; i++)
            {
                DownIndexList[i].AddRange(UpIndexList[i]);
                DownIndexList[i].Reverse();
            }

            multiPatchGeometryCollection[TINID] = new MultiPatchClass();
            IPointCollection[] ringPointCollection = new IPointCollection[4];
            int[] BoundryCount = new int[4];
            for (int i = 0; i < 4; i++)
            {
                BoundryCount[i] = UpIndexList[i].Count;
                ringPointCollection[i] = new RingClass();
            }

            for (int i = 0; i < 4; i++)
            {
                //up
                for (int j = 0; j < BoundryCount[i]; j++)
                {
                    IPoint Upoint = new PointClass();
                    Material.MakeZAware(Upoint);
                    orderNodeU[(int)(UpIndexList[i])[j]].QueryAsPoint(Upoint);
                    ringPointCollection[i].AddPoint(Upoint, ref _missing, ref _missing);
                }
                //down
                for (int j = 0; j < BoundryCount[i]; j++)
                {
                    IPoint Dpoint = new PointClass();
                    Material.MakeZAware(Dpoint);
                    orderNodeD[(int)(DownIndexList[i])[j]].QueryAsPoint(Dpoint);
                    ringPointCollection[i].AddPoint(Dpoint, ref _missing, ref _missing);
                }
                //最后一点
                IPoint LastPoint = new PointClass();
                Material.MakeZAware(LastPoint);
                orderNodeU[(int)(UpIndexList[i])[0]].QueryAsPoint(LastPoint);
                ringPointCollection[i].AddPoint(LastPoint, ref _missing, ref _missing);
            }
            for (int j = 0; j < 4; j++)
            {
                multiPatchGeometryCollection[TINID].AddGeometry(ringPointCollection[j] as IGeometry, ref _missing, ref _missing);
            }

            ////添加上下曲面
            //up
            int triangleCountU = tinU.TriangleCount;
            for (int i = 0; i < triangleCountU; i++)
            {
                bool add = true;
                ITinTriangle triangleU = tinU.GetTriangle(i + 1);
                if (triangleU.IsEmpty)
                    continue;
                IPointCollection trianglePointCollection = new TrianglesClass(); ;
                for (int j = 0; j < 3; j++)
                {
                    IPoint tp = new PointClass();
                    Material.MakeZAware(tp);
                    (triangleU.get_Node(j)).QueryAsPoint(tp);
                    if (tp.Z.Equals(Double.NaN))
                    {
                        add = false;
                        break;
                    }
                    trianglePointCollection.AddPoint(tp, ref _missing, ref _missing);
                }
                if (add == true)
                {
                    multiPatchGeometryCollection[TINID].AddGeometry(trianglePointCollection as IGeometry, ref _missing, ref _missing);
                }
            }

            //down
            int triangleCountD = tinD.TriangleCount;
            for (int i = 0; i < triangleCountD; i++)
            {
                bool add = true;
                ITinTriangle triangleD = tinD.GetTriangle(i + 1);
                if (triangleD.IsEmpty)
                    continue;
                IPointCollection trianglePointCollection = new TrianglesClass(); ;
                for (int j = 0; j < 3; j++)
                {
                    IPoint tp = new PointClass();
                    Material.MakeZAware(tp);
                    (triangleD.get_Node(j)).QueryAsPoint(tp);
                    if (tp.Z.Equals(Double.NaN))
                    {
                        add = false;
                        break;
                    }
                    trianglePointCollection.AddPoint(tp, ref _missing, ref _missing);
                }
                if (add == true)
                {
                    multiPatchGeometryCollection[TINID].AddGeometry(trianglePointCollection as IGeometry, ref _missing, ref _missing);
                }
            }
        }

        public void CreateFeature()
        {
            //移动放缩矢量箭头
            IVector3D VerticalVector = Material.ConstructVector3D(0, 0, 1);
            IPoint originPoint = Material.ConstructPoint3D(0, 0, 0);
            Material.MakeZAware(originPoint as IGeometry);

            //创建所有要素
            for (int i = 0; i < VectorModel.TinCount; i++)
            {
                IFeature BarycentreFeature = BarycentreFeatureClass.CreateFeature();
                BarycentreFeature.Shape = BARYCENTRE[i];
                int index1 = BarycentreFeature.Fields.FindField("X");
                BarycentreFeature.set_Value(index1, BARYCENTRE[i].X);
                int index2 = BarycentreFeature.Fields.FindField("Y");
                BarycentreFeature.set_Value(index2, BARYCENTRE[i].Y);
                int index3 = BarycentreFeature.Fields.FindField("Z");
                BarycentreFeature.set_Value(index3, BARYCENTRE[i].Z);
                int index4 = BarycentreFeature.Fields.FindField("Probability");
                BarycentreFeature.set_Value(index4, 0.95 - (0.9 / 49) * i);
                BarycentreFeature.Store();

                if (i == 0)
                {
                    Ctxt.WriteLine("GOCAD VSet" + "\r\n"
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
                }
                Ctxt.Write("PVRTX " + i + " ");
                Ctxt.Write(BARYCENTRE[i].X + " ");
                Ctxt.Write(BARYCENTRE[i].Y + " ");
                Ctxt.Write(BARYCENTRE[i].Z + " ");
                if (i == VectorModel.TinCount - 1)
                {
                    Ctxt.WriteLine("0.5" + " ");
                }
                else
                {
                    Ctxt.WriteLine(0.95 - (0.9 / 49) * i + " ");
                }

                if (i == VectorModel.TinCount - 1)
                {
                    Ctxt.Write("END");
                    Ctxt.Close();
                }


                ////创建重力线
                IFeature VectorFeature = GravityVectorFeatureClass.CreateFeature();
                IPoint endPoint = new PointClass();
                Material.MakeZAware(endPoint as IGeometry);
                //根据体积设置末端点
                endPoint.X = BARYCENTRE[i].X;
                endPoint.Y = BARYCENTRE[i].Y;
                endPoint.Z = BARYCENTRE[i].Z - VOLUME[i] * gLength;                                        //
                //添加两点构成向量
                IPointCollection VectorPointCollection = new PolylineClass();
                Material.MakeZAware(VectorPointCollection as IGeometry);
                VectorPointCollection.AddPoint(BARYCENTRE[i], ref _missing, ref _missing);
                VectorPointCollection.AddPoint(endPoint, ref _missing, ref _missing);
                VectorFeature.Shape = VectorPointCollection as IGeometry;
                int index5 = VectorFeature.Fields.FindField("Gravity");
                VectorFeature.set_Value(index5, VOLUME[i]);
                int index6 = VectorFeature.Fields.FindField("Probability");
                VectorFeature.set_Value(index6, 0.95 - (0.9 / 49) * i);
                VectorFeature.Store();

                ////添加矢量箭头
                //得到重力矢量
                IVector3D Gvector = Material.CreateVector3DTwoPoints(endPoint, BARYCENTRE[i]);
                IGeometry Arrow = Material.GetArrow();
                double Inclination = Gvector.Inclination;
                //计算与竖向矢量夹角
                double degreesOfRotation = Math.Acos((Gvector.DotProduct(VerticalVector)) / ((Gvector.Magnitude) * (VerticalVector.Magnitude)));
                ITransform3D transform3D = Arrow as ITransform3D;
                //放缩
                transform3D.Scale3D(originPoint, XScale, YScale, ZScale);
                //转动
                if (degreesOfRotation != 0)
                {
                    double angleOfRotationInRadians = degreesOfRotation;
                    IVector3D axisOfRotationVector3D = new Vector3DClass();
                    axisOfRotationVector3D.XComponent = 1;
                    axisOfRotationVector3D.YComponent = 0;
                    axisOfRotationVector3D.ZComponent = 0;
                    transform3D.RotateVector3D(axisOfRotationVector3D, angleOfRotationInRadians);
                }
                //平移
                if (endPoint.IsEmpty)
                    continue;
                transform3D.Move3D(endPoint.X - originPoint.X, endPoint.Y - originPoint.Y, endPoint.Z - originPoint.Z);
                IFeature ArrowFeature = ArrowFeatureClass.CreateFeature();
                ArrowFeature.Shape = Arrow as IMultiPatch;
                int index7 = ArrowFeature.Fields.FindField("Probability");
                ArrowFeature.set_Value(index7, 0.95 - (0.9 / 49) * i);
                ArrowFeature.Store();
                //墙壁
                IFeature EnCloseFeature = EncloseFeatureClass[i].CreateFeature();
                EnCloseFeature.Shape = multiPatchGeometryCollection[i] as IMultiPatch;
                int index17 = EnCloseFeature.Fields.FindField("Probability");
                EnCloseFeature.set_Value(index7, 0.95 - (0.9 / 49) * i);
                EnCloseFeature.Store();

                Vtxt.WriteLine(VOLUME[i]);
                if (i == VectorModel.TinCount - 1)
                {
                    Vtxt.Close();
                }
            }
        }

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
                //Create a user-defined double field for "Probability"
                IFieldsEdit fieldsEdit148 = (IFieldsEdit)fields;
                IField incomeField148 = new FieldClass();
                IFieldEdit incomeFieldEdit148 = (IFieldEdit)incomeField148;
                incomeFieldEdit148.AliasName_2 = "Probability";
                incomeFieldEdit148.Editable_2 = true;
                incomeFieldEdit148.IsNullable_2 = false;
                incomeFieldEdit148.Name_2 = "Probability";
                incomeFieldEdit148.Precision_2 = 2;
                incomeFieldEdit148.Scale_2 = 5;
                incomeFieldEdit148.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit148.AddField(incomeField148);
            }
            //If it is a polyline, add Field
            if (geometryType == esriGeometryType.esriGeometryPolyline)
            {
                //Create a user-defined double field for "Magnitude"
                IFieldsEdit fieldsEdit = (IFieldsEdit)fields;
                IField incomeField = new FieldClass();
                IFieldEdit incomeFieldEdit = (IFieldEdit)incomeField;
                incomeFieldEdit.AliasName_2 = "Gravity";
                incomeFieldEdit.Editable_2 = true;
                incomeFieldEdit.IsNullable_2 = false;
                incomeFieldEdit.Name_2 = "Gravity";
                incomeFieldEdit.Precision_2 = 2;
                incomeFieldEdit.Scale_2 = 5;
                incomeFieldEdit.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit.AddField(incomeField);
                //Create a user-defined double field for "Probability"
                IFieldsEdit fieldsEdit148 = (IFieldsEdit)fields;
                IField incomeField148 = new FieldClass();
                IFieldEdit incomeFieldEdit148 = (IFieldEdit)incomeField148;
                incomeFieldEdit148.AliasName_2 = "Probability";
                incomeFieldEdit148.Editable_2 = true;
                incomeFieldEdit148.IsNullable_2 = false;
                incomeFieldEdit148.Name_2 = "Probability";
                incomeFieldEdit148.Precision_2 = 2;
                incomeFieldEdit148.Scale_2 = 5;
                incomeFieldEdit148.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit148.AddField(incomeField148);
            }

            if (geometryType == esriGeometryType.esriGeometryMultiPatch)
            {
                //Create a user-defined double field for "Probability"
                IFieldsEdit fieldsEdit148 = (IFieldsEdit)fields;
                IField incomeField148 = new FieldClass();
                IFieldEdit incomeFieldEdit148 = (IFieldEdit)incomeField148;
                incomeFieldEdit148.AliasName_2 = "Probability";
                incomeFieldEdit148.Editable_2 = true;
                incomeFieldEdit148.IsNullable_2 = false;
                incomeFieldEdit148.Name_2 = "Probability";
                incomeFieldEdit148.Precision_2 = 2;
                incomeFieldEdit148.Scale_2 = 5;
                incomeFieldEdit148.Type_2 = esriFieldType.esriFieldTypeDouble;
                fieldsEdit148.AddField(incomeField148);
            }
            return fields;
        }

        public void CreateFeatureClass()
        {
            string Barycentre = "Barycentre";
            string GravityVector = "GravityVector";
            string Arrow = "Arrow";
            string[] Enclose = new string[VectorModel.TinCount];
            Ctxt = new StreamWriter(Cpath);
            Vtxt = new StreamWriter(Vpath);


            //传递坐标系
            IGeoDataset geoDataset = tinD as IGeoDataset;
            ISpatialReference spatialReference = geoDataset.SpatialReference;

            //创建要素类
            BarycentreFields = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPoint);
            BarycentreFeatureClass = CreateStandaloneFeatureClass(Workspace, Barycentre, BarycentreFields);

            GravityVectorFields = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryPolyline);
            GravityVectorFeatureClass = CreateStandaloneFeatureClass(Workspace, GravityVector, GravityVectorFields);

            ArrowFields = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryMultiPatch);
            ArrowFeatureClass = CreateStandaloneFeatureClass(Workspace, Arrow, ArrowFields);

            EncloseFields = CreateFieldsCollectionForFeatueClass(spatialReference, esriGeometryType.esriGeometryMultiPatch);
            for (int i = 0; i < VectorModel.TinCount; i++)
            {
                Enclose[i] = "Enclose" + i;
                EncloseFeatureClass[i] = CreateStandaloneFeatureClass(Workspace, Enclose[i], EncloseFields);
            }

        }
        #endregion Prepare

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
        #endregion
    }
}
