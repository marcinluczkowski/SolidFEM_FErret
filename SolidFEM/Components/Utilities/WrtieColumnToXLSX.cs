using System;
using System.Collections.Generic;
using ClosedXML.Excel;
using SolidFEM.Classes;
using Grasshopper.Kernel;

namespace SolidFEM.Utilities
{
    public class WriteColumnXLSX : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public WriteColumnXLSX()
            : base("WriteColumnXLSX", "col2xlsx",
                "Preview nodal stresses and displacements. ",
                "SmartMesh", "DataHandling")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Path", "path", "Path of Excel file", GH_ParamAccess.item); // 0
            pManager.AddTextParameter("Sheet", "sheet", "Sheet to write to. If empty, The first sheet is used",
                GH_ParamAccess.item); //1
            pManager.AddIntegerParameter("Column", "col", "Position of column to write data to", GH_ParamAccess.item); // 2
            pManager.AddGenericParameter("List", "lst", "List of data to be written to column", GH_ParamAccess.list); // 3

            pManager[1].Optional = true; // do not have to specify sheet
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            // No output from this file
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // -- input -- 
            string path = ""; // 0
            string sheet = "Sheet1"; // 1
            int col = 0;
            List<double> colLst = new List<double>(); // 3

            if (!DA.GetData(0, ref path)) return; // 0
            DA.GetData(1, ref path); // 1
            DA.GetData(2, ref col);
            if (!DA.GetDataList(2, colLst)) return; // 3


            // -- solve --

            // try to open the file
            XLWorkbook wb = new XLWorkbook();
            try
            {
                wb = new XLWorkbook(path);
            }
            catch (Exception e)
            {
                Console.WriteLine(e);
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, e.ToString());
            }

            // Specify sheet to write to
            IXLWorksheet ws;
            if (!wb.TryGetWorksheet(sheet, out ws))
            {
                ws = wb.Worksheet(0);
            }

            // Write column to sheet. Start on index 1 to preserve column header. 
            for (int i = 1; i < colLst.Count + 1; i++)
            {
                ws.Column(col).Cell(i).SetValue(colLst[i]);
            }

            // -- output
            wb.SaveAs(path); // save the modified file. 
            AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "Column successfully written to:" + path + "in Sheet: " + sheet);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("d4554d2f-2b2f-4a35-9bf2-a0c017d17b82"); }
        }
    }
}