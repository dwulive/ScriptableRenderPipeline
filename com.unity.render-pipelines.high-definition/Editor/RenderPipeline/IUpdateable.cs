using UnityEditor.AnimatedValues;
using UnityEngine.Events;

namespace UnityEditor.Rendering.HDPipeline
{
    public interface IUpdateable<T>
    {
        void Update(T v);
    }
}
