# -*- coding: utf-8 -*-

from django.shortcuts import render
from django.views.decorators import csrf
from pc import Daoxian


# 接收POST请求数据
def search_post(request):
    ctx = {}
    if request.POST:
        ctx['result'] = request.POST
        try:
            dxw = Daoxian(request.POST['info1'], request.POST['info2'], request.POST['info3'])
            dxw.check_info()
            dxw.init_params()
            dxw.init_weight()
            dxw.cal()
            ctx['result'] = dxw.params_to_html_th()
        except:
            ctx['result'] = 'error input'
    return render(request, "dxw_solve.html", ctx)