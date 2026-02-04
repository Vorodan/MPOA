using System;

namespace MPOA
{
    public class AlgParameter
    {
        public String Name;
        public Object Value;

        public AlgParameter(String name, Object value) { this.Name = name; this.Value = value; }
        public AlgParameter() { }
    }
}
