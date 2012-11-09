#include <memory>
#include <sstream>
#include "network_package.h"
#include "network_wrapper.h"
#include <cassert>

using namespace std;

//Static Memeber Initialization
int NetPackage::last_group_id = 0;

//Static Methods
int NetPackage::get_new_group_id()
{
	return ++last_group_id;
}

//Non static methods
NetPackage::NetPackage(PackType _type, int _group, int _run_id, const string &desc_str)
	: type(_type), group(_group), run_id(_run_id)
{
	memset(desc, '\0', DESC_LEN);
	strncpy(desc, desc_str.c_str(), DESC_LEN-1);
	data_len = 1;
}
void NetPackage::reset(PackType _type, int _group, int _run_id, const string &_desc)
{
	type = _type;
	group = _group;
	run_id = _run_id;
	memset(desc, '\0', DESC_LEN);
	strncpy(desc, _desc.c_str(), DESC_LEN-1);
	data.clear();
}

int NetPackage::send(int sockfd, const void *data, unsigned long data_len_l)
{
	int n;
	unsigned long buf_sz = 0;
	//calculate the size of buffer
	buf_sz = sizeof(buf_sz);
	buf_sz += sizeof(type);
	buf_sz += sizeof(group);
	buf_sz += sizeof(run_id);
	buf_sz += sizeof(desc);
	buf_sz += data_len_l;
	//pack information into buffer
	//unique_ptr<char[]> buf(new char[buf_sz]);
	vector<char> buf;
	buf.resize(buf_sz, '\0');
	size_t i_start = 0;
	w_memcpy_s(&buf[i_start], buf_sz-i_start, &buf_sz, sizeof(buf_sz));
	i_start += sizeof(buf_sz);
	w_memcpy_s(&buf[i_start], buf_sz-i_start, &type, sizeof(type));
	i_start += sizeof(type);
	w_memcpy_s(&buf[i_start], buf_sz-i_start, &group, sizeof(group));
	i_start += sizeof(group);
	w_memcpy_s(&buf[i_start], buf_sz-i_start, &run_id, sizeof(run_id));
	i_start += sizeof(run_id);
	w_memcpy_s(&buf[i_start], buf_sz-i_start, desc, sizeof(desc));
	i_start += sizeof(desc);
	if (data_len_l > 0) {
		w_memcpy_s(&buf[i_start], buf_sz-i_start, data, data_len_l);
		i_start += data_len_l;
	}
	if (i_start!=buf_sz) {
		cerr << "NetPackage::send error: could ony send" << i_start
			<< " out of " << buf_sz << "bytes" << endl;
	}
	assert (i_start==buf_sz);
	n = w_sendall(sockfd, buf.data(), &buf_sz);
	return n;  // return 0 on sucess or -1 on failure
}

int  NetPackage::recv(int sockfd)
{
	int n;
	unsigned long header_sz = 0;
	unsigned long buf_sz = 0;
	size_t i_start = 0;
	//get header (ie size, seq_id, id and name)
	header_sz = sizeof(buf_sz) + sizeof(type) + sizeof(group) + sizeof(run_id) + sizeof(desc);
	vector<char> header_buf;
	header_buf.resize(header_sz, '\0');
	n = w_recvall(sockfd, &header_buf[0], &header_sz);
	if(n>0) {
		assert(header_sz==header_buf.size());
		i_start = 0;
		w_memcpy_s(&buf_sz, sizeof(buf_sz), &header_buf[i_start], sizeof(buf_sz));
		i_start += sizeof(buf_sz);
		w_memcpy_s(&type, sizeof(type), &header_buf[i_start], sizeof(type));
		i_start += sizeof(type);
		w_memcpy_s(&group, sizeof(group), &header_buf[i_start], sizeof(group));
		i_start += sizeof(group);
		w_memcpy_s(&run_id, sizeof(run_id), &header_buf[i_start], sizeof(run_id));
		i_start += sizeof(run_id);
		w_memcpy_s(&desc, sizeof(desc), &header_buf[i_start], sizeof(desc));
		i_start += sizeof(desc);
		desc[DESC_LEN-1] = '\0';
		//get data
		data_len = buf_sz - i_start;
		data.resize(data_len, '\0');
		if (data_len > 0) {
			n = w_recvall(sockfd, &data[0], &data_len);
			assert(data_len==buf_sz-header_sz);
		}
	}
	if (n> 1) {n=1;}
	return n;  // -1 on failure, 0 on a close connection or 1 on success
}

void NetPackage::print_header(std::ostream &fout)
{
	fout << "NetPackage: type = " << int(type) <<", group = " << group << ", run_id = " << run_id << ", description = " << desc << 
		", data package size = " << data.size() << endl; 
}
