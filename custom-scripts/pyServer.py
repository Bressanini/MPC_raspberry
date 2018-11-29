import SimpleHTTPServer
import SocketServer
import os.path, sys
import datetime
import sys
import BaseHTTPServer
import datetime
import psutil
import time
import os
import platform
import subprocess

from SimpleHTTPServer import SimpleHTTPRequestHandler

#Formating string. The %s will be substituted by system data.
html = "<!DOCTYPE html>  \
		<html> \
			<body> \
				<h1>System Information</h1> \
				<p>Local Data: %s</p> "
			

#Handle de http requests.
class MyRequestHandler(SimpleHTTPServer.SimpleHTTPRequestHandler):
	#Build the html file for every connection.
	def makeindex(self):
		html = "<!DOCTYPE html>  \
		<html> \
            <head>\
                <meta charset=\"UTF-8\">\
            </head>\
            <style>\
                body{background-color: #000000}\
                table{background-color: #c9c9c9}\
                table, th, td {\
                    border: 1px solid black;\
                }\
                th, td{background-color: #cdd2d8}\
                #text{color:WHITE}\
            </style>\
			<body><center> \
				<h1 id=\"text\">Informacoes do Sistema</h1>\
                <p id=\"text\"><b>Grupo: </b> Leticia Marchi, Gabriel Bressanini</p>"
		with open('/proc/uptime', 'r') as f:
			uptime_seconds = float(f.readline().split()[0])
		with open('/proc/cpuinfo') as f:
				       for line in f:           
				               if line.strip():
				                       if line.rstrip('\n').startswith('model name'):
				                               model_name = line.rstrip('\n').split(':')[1]
		dateTime = datetime.datetime.now().time()
		timea = time.strftime("%c")
		processes = filter(lambda p: psutil.Process(p).name() != "python", psutil.pids())
		cpu_informacao = model_name
		cpu_vel = psutil.cpu_freq();
		cpu_perc = psutil.cpu_percent(0);
		mem = psutil.virtual_memory()
		mem_total = mem.total /(1024*1024);
		mem_used = mem.used /(1024*1024);
		os_sys = platform.system();
		os_rel = platform.release();
		os_ver = platform.version();
		html = str(html)+ "<table>\
				        <tr>\
				             <td><p><b>Data e hora do sistema:</b></td> <td align=\"center\">{code}</p></td>\
				        </tr>\
				        <tr>\
				             <td><p><b>Uptime:</b></td> <td align=\"center\">{code1} segundos</p></td>\
				        </tr>\
				        <tr>\
				             <td><p><b>Informacoes dobre a CPU:</b></td> <td align=\"center\">{code2}</p></td>\
				        </tr>\
				        <tr>\
				             <td><p><b>Capacidade ocupada do Processador:</b></td> <td align=\"center\">{code3}%</p></td>\
				        </tr>\
				        <tr>\
				             <td><p><b>Quantidade de memoria RAM total:</b></td> <td align=\"center\">{code4}MB</p></td>\
				        </tr>\
				        <tr>\
				             <td><p><b>Quantidade de memoria RAM utilizada:</b></td> <td align=\"center\">{code5}MB</td>\
				        </tr>\
				        <tr>\
				             <td><p><b>Versao do Sistema:</b></td> <td align=\"center\"><p>{code6}</p><p>{code7}</p><p>{code8}</p></p></td>\
				        </tr>\
				        </table>\
				             <p id=\"text\"><b>Lista de processos em execucao:</b></p>".format(code=timea,code1=uptime_seconds,code2=cpu_informacao,code3=cpu_perc,code4=mem_total,code5=mem_used,code6=os_sys,code7=os_rel,code8=os_ver)
		html = str(html) + "<table>\
		  <tr>\
		    <th>PID</th>\
		    <th>Nome</th>"
		for proc in psutil.process_iter():
			name = proc.name()
			pid1 = proc.pid
			html = str(html) + "</tr>\
				    <td>{code2}</td>\
				    <td align=\"center\">{code}</td>\
				     </tr>".format(code=name,code2=pid1)

		html = str(html) + "</table></center></body></html>"
		f = open("index.html", "w")
		f.write(html)
		f.close()
		return
	
	#Method http GET.	
	def do_GET(self):
		self.makeindex()
		self.path = '/index.html'
		return SimpleHTTPServer.SimpleHTTPRequestHandler.do_GET(self)

#Start the webserver.
Handler = MyRequestHandler
server = SocketServer.TCPServer(('0.0.0.0', 8080), Handler)
server.serve_forever()	
