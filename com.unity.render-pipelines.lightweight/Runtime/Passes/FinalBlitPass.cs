namespace UnityEngine.Rendering.LWRP
{
    /// <summary>
    /// Copy the given color target to the current camera target
    ///
    /// You can use this pass to copy the result of rendering to
    /// the camera target. The pass takes the screen viewport into
    /// consideration.
    /// </summary>
    internal class FinalBlitPass : ScriptableRenderPass
    {
        const string m_ProfilerTag = "Final Blit Pass";
        RenderTargetHandle m_Source;
        Material m_BlitMaterial;
        TextureDimension m_TargetDimension;
        bool m_ClearBlitTarget;
        bool m_IsMobileOrSwitch;
        Rect m_PixelRect;

        public FinalBlitPass(RenderPassEvent evt, Material blitMaterial)
        {
            m_BlitMaterial = blitMaterial;
            renderPassEvent = evt;
        }

        /// <summary>
        /// Configure the pass
        /// </summary>
        /// <param name="baseDescriptor"></param>
        /// <param name="colorHandle"></param>
        /// <param name="clearBlitTarget"></param>
        /// <param name="pixelRect"></param>
        public void Setup(RenderTextureDescriptor baseDescriptor, RenderTargetHandle colorHandle, bool clearBlitTarget = false, Rect pixelRect = new Rect())
        {
            m_Source = colorHandle;
            m_TargetDimension = baseDescriptor.dimension;
            m_ClearBlitTarget = clearBlitTarget;
            m_IsMobileOrSwitch = Application.isMobilePlatform || Application.platform == RuntimePlatform.Switch;
            m_PixelRect = pixelRect;
        }

        /// <inheritdoc/>
        public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
        {
            if (m_BlitMaterial == null)
            {
                Debug.LogErrorFormat("Missing {0}. {1} render pass will not execute. Check for missing reference in the renderer resources.", m_BlitMaterial, GetType().Name);
                return;
            }

            bool requiresSRGBConvertion = Display.main.requiresSrgbBlitToBackbuffer;
            bool killAlpha = renderingData.killAlphaInFinalBlit;

            CommandBuffer cmd = CommandBufferPool.Get(m_ProfilerTag);

            if (requiresSRGBConvertion)
                cmd.EnableShaderKeyword(ShaderKeywordStrings.LinearToSRGBConversion);
            else
                cmd.DisableShaderKeyword(ShaderKeywordStrings.LinearToSRGBConversion);

            if (killAlpha)
                cmd.EnableShaderKeyword(ShaderKeywordStrings.KillAlpha);
            else
                cmd.DisableShaderKeyword(ShaderKeywordStrings.KillAlpha);

            ref CameraData cameraData = ref renderingData.cameraData
            ;
            // Use default blit for XR as we are not sure the UniversalRP blit handles stereo.
            // The blit will be reworked for stereo along the XRSDK work.
            var pipe = GraphicsSettings.renderPipelineAsset as LightweightRenderPipelineAsset;

            Material blitMaterial = (pipe.blitMaterial != null) ? pipe.blitMaterial : (cameraData.isStereoEnabled) ? null : m_BlitMaterial;
            cmd.SetGlobalTexture("_BlitTex", m_Source.Identifier());
            if (cameraData.isStereoEnabled || cameraData.isSceneViewCamera || cameraData.isDefaultViewport)
            {
                // This set render target is necessary so we change the LOAD state to DontCare.
                cmd.SetRenderTarget(BuiltinRenderTextureType.CameraTarget,
                    RenderBufferLoadAction.DontCare,
				   cameraData.isSceneViewCamera ? RenderBufferStoreAction.Store : RenderBufferStoreAction.DontCare,
                    RenderBufferLoadAction.DontCare,
                    cameraData.isSceneViewCamera ? RenderBufferStoreAction.Store : RenderBufferStoreAction.DontCare);

                // Clearing render target is cost free on mobile and it avoid tile loading
                // UGEN
                // Must clear for Screen Spawn canvases.  Don't need to clear color, but I this is free on mobile tiled renderers.
                if (true ) // m_IsMobileOrSwitch)
                    cmd.ClearRenderTarget(true, true, Color.black);
				
                cmd.Blit(m_Source.Identifier(), BuiltinRenderTextureType.CameraTarget, blitMaterial);
            }
            else
            {
                // TODO: Final blit pass should always blit to backbuffer. The first time we do we don't need to Load contents to tile.
                // We need to keep in the pipeline of first render pass to each render target to propertly set load/store actions.
                // meanwhile we set to load so split screen case works.
                SetRenderTarget(
                    cmd,
                    BuiltinRenderTextureType.CameraTarget,
                    m_ClearBlitTarget ? RenderBufferLoadAction.DontCare : RenderBufferLoadAction.Load,
                    RenderBufferStoreAction.DontCare,
                    RenderBufferLoadAction.DontCare,
                    RenderBufferStoreAction.DontCare,
                    (m_ClearBlitTarget ? ClearFlag.Color : ClearFlag.None)|ClearFlag.Depth,
                    Color.black,
                    m_TargetDimension);

                Camera camera = cameraData.camera;
                cmd.SetViewProjectionMatrices(Matrix4x4.identity, Matrix4x4.identity);
                cmd.SetViewport(m_PixelRect != Rect.zero ? m_PixelRect : camera.pixelRect);
                cmd.DrawMesh(RenderingUtils.fullscreenMesh, Matrix4x4.identity, blitMaterial);
    
            }
            if (pipe.uiCamera != null)
            {
                var worldToCameraMatrix = cameraData.camera.worldToCameraMatrix;
                var projectionMatrix = pipe.uiCamera.projectionMatrix;
                cmd.SetViewProjectionMatrices(worldToCameraMatrix, projectionMatrix);
                LightweightRenderPipeline.ModifyCameraShaderConstants(ref worldToCameraMatrix, ref projectionMatrix, cameraData.camera.transform.position,cmd);
            }
            CoreUtils.SetKeyword(cmd, ShaderKeywordStrings.MainLightShadows, false);

            context.ExecuteCommandBuffer(cmd);
            CommandBufferPool.Release(cmd);
        }
    }
}
